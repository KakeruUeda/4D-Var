#include "DataGrid.h"

void DataGrid::weighted_average_time(const int iv, const int t)
{
  if(cfd_time_steps_id[t].size() == 0) return;

  auto get_node_values = [&](const int it, const int t) {
    mt1d.xCurrent(0) = (it - 1) * dt_cfd;
    mt1d.xCurrent(1) = it * dt_cfd;

    for(int d = 0; d < 3; d++) {
      vel_current_time(0, d) = voxel(iv).v_cfd_refined(it - 1, d);
      vel_current_time(1, d) = voxel(iv).v_cfd_refined(it, d);
    }
  };

  auto evaluate = [&](double &s_gp, vector<double> &v_gp, const int iv, const int t) {
    weight_integral_time += s_gp * mt1d.vol;
    for(int d = 0; d < 3; d++) {
      voxel(iv).v_cfd(t, d) += s_gp * v_gp[d] * mt1d.vol;
    }
  };

  auto average = [&](const int iv, const int t) {
    for(int i1 = 0; i1 < 2; i1++) {
      mt1d.setShapesInGauss(gauss, i1);
      mt1d.setFactorsInGauss(gauss, i1, dt_cfd);
      auto s_gp = mt1d.getScalarValueGP(smoothing_time);
      auto v_gp = mt1d.getVectorValueGP(vel_current_time);
      evaluate(s_gp, v_gp, iv, t);
    }
  };

  auto compute_smoothing_time = [&](const int t) {
    double sigma = 0.002;

    auto kai = [](double value, double center, double d, double a) {
      return 1 / (1 + exp(-(value - (center - (d / 2e0))) / a)) - 1 / (1 + exp(-(value - (center + (d / 2e0))) / a));
    };

    auto compute_kai = [&](int p) { return kai(mt1d.xCurrent(p), t * dt_mri, dt_mri, sigma); };

    smoothing_time(0) = compute_kai(0);
    smoothing_time(1) = compute_kai(1);
  };

  weight_integral_time = 0e0;

  int it;

  for(int i = 0; i < cfd_time_steps_id[t].size() + 1; i++) {
    if(i == cfd_time_steps_id[t].size()) {
      it = cfd_time_steps_id[t][i - 1] + 1;
    } else {
      it = cfd_time_steps_id[t][i];
    }

    if(it - 1 < 0) continue;
    if(it > n_cfd_step - 1) continue;

    get_node_values(it, t);
    compute_smoothing_time(t);
    average(iv, t);
  }

  if(weight_integral_time == 0) {
    if(mpi.myId == 0) std::cout << "weight_integral_time == 0. Exit." << std::endl;
    MPI_Finalize();
    exit(1);
  }

  for(int d = 0; d < 3; d++) {
    voxel(iv).v_cfd(t, d) /= weight_integral_time;
  }
}

void DataGrid::spline_interpolate_time(const int iv, const int itm)
{
  std::vector<double> x, y1, y2, y3;

  for(int itc = 0; itc < n_cfd_step; itc++) {
    double p = itc * dt_cfd;
    x.push_back(p);
    y1.push_back(voxel(iv).v_cfd_refined(itc, 0));
    y2.push_back(voxel(iv).v_cfd_refined(itc, 1));
    y3.push_back(voxel(iv).v_cfd_refined(itc, 2));
  }

  vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
  vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
  vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);

  double p = itm * dt_mri;
  voxel(iv).v_cfd(itm, 0) = Spline::evaluate(cf_x, p);
  voxel(iv).v_cfd(itm, 1) = Spline::evaluate(cf_y, p);
  voxel(iv).v_cfd(itm, 2) = Spline::evaluate(cf_z, p);
}
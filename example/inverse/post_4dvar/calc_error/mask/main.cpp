#include <iostream>
#include <fstream>
#include <vector>
#include <string>

static const int DX[6] = {1, -1,  0,  0,  0,  0};
static const int DY[6] = {0,  0,  1, -1,  0,  0};
static const int DZ[6] = {0,  0,  0,  0,  1, -1};

inline int idx3D(int i, int j, int k, int nx, int ny, int nz)
{
    return i + nx * (j + ny * k);
}

int main()
{
    const int nx = 40;
    const int ny = 40;
    const int nz = 27;
    const int N  = nx * ny * nz;

    // 1. ファイルから読み込む (入力配列: mask_in)
    std::vector<double> mask_in(N, 0.0);
    {
        std::string filename = "mask.dat";
        std::ifstream ifs(filename.c_str());
        if(!ifs){
            std::cerr << "Error: Could not open " << filename << std::endl;
            return -1;
        }
        for(int i = 0; i < N; i++){
            if(!(ifs >> mask_in[i])){
                std::cerr << "Error: Not enough data in " << filename << std::endl;
                return -1;
            }
        }
        ifs.close();
    }

    // 3. まず「界面 (0 < val < 1)」を0にする (mask_inを見て判断し、mask_outに書き込む)
    for(int k = 0; k < nz; k++){
        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
                int idx = idx3D(i, j, k, nx, ny, nz);
                double val = mask_in[idx];
                if(val > 0.0 && val < 1.0){
                    mask_in[idx] = 0.0; 
                }
            }
        }
    }

    // 2. 出力用の配列 (mask_out) を最初は mask_in と同じにしておく
    std::vector<double> mask_out = mask_in;

    // 4. 一回だけ「0 に隣接する流体(=1)」を0にする
    //    ※ "mask_in" を参照して近傍に0があるかを見て、"mask_out" に書き込み
    //    ※ 近傍に「既に界面のせいで0になった」セルがあるかを確かめたい場合、
    //       ここでは mask_in をそのまま使うか、(3) で書き換えた mask_out を使うかで挙動が違う
    //       「界面が0になったあと、そこに隣接する流体を消す」なら (3) で書き換えた mask_out を見るのが自然
    //       ただし「元々の mask_in が 0 だった(固体)」も考慮するなら mask_in も確認が必要
    //    →ここは「界面を0にした後の状態」を見たいはずなので、mask_out を基準にする例を紹介
    {
        // mask_out をもとに「どこが 0 になっているか」判定し、その近傍にいた mask_in=1 の場所を mask_out=0 にする
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++){
                for(int i = 0; i < nx; i++){
                    int idx = idx3D(i, j, k, nx, ny, nz);
                    // 元々流体(=1) だったか？ (mask_in を見てもいいし、mask_out を見てもいいが、ここでは「最初の状態」を見たいなら mask_in)
                    if(mask_in[idx] == 1.0) {
                        // 近傍に0があるか？ (0はmask_outでチェック)
                        bool hasZeroNeighbor = false;
                        for(int n = 0; n < 6; n++){
                            int ii = i + DX[n];
                            int jj = j + DY[n];
                            int kk = k + DZ[n];
                            if(ii < 0 || ii >= nx) continue;
                            if(jj < 0 || jj >= ny) continue;
                            if(kk < 0 || kk >= nz) continue;
                            
                            int nIdx = idx3D(ii, jj, kk, nx, ny, nz);
                            if(mask_in[nIdx] == 0.0) {
                                hasZeroNeighbor = true;
                                break;
                            }
                        }
                        if(hasZeroNeighbor){
                            // 一回だけの処理なので mask_out を 0 に更新
                            mask_out[idx] = 0.0;
                        }
                    }
                }
            }
        }
    }

    // 5. 結果をファイルに書き出したい場合 (例)
    {
        std::ofstream ofs("mask_x.dat");
        if(!ofs){
            std::cerr << "Error: Could not open output file." << std::endl;
            return -1;
        }
        for(int i = 0; i < N; i++){
            ofs << mask_out[i] << "\n";
        }
        ofs.close();
    }

    std::cout << "Done!" << std::endl;
    return 0;
}

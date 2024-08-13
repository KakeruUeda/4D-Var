/**
 * @file Config.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderInterface::readBasicInfo()
{
  std::string str, base_label, label;
  base_label = "/Base";
  label = base_label + "/dim";

  if(!ptr->tp.getInspectedValue(label, ptr->dim)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/outputDir";
  if(!ptr->tp.getInspectedValue(label, ptr->outputDir)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderInterface::readPhysicalInfo()
{
  std::string str, base_label, label;

  base_label = "/PysicalParameter";
  label = base_label + "/rho";
  if(!ptr->tp.getInspectedValue(label, ptr->rho)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/mu";
  if(!ptr->tp.getInspectedValue(label, ptr->mu)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderInterface::readTimeInfo()
{
  std::string str, base_label, label;

  base_label = "/TimeParameter";
  label = base_label + "/dt";
  if(!ptr->tp.getInspectedValue(label, ptr->dt)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/timeMax";
  if(!ptr->tp.getInspectedValue(label, ptr->timeMax)) {
    throw std::runtime_error(label + " is not set");
  }

  std::string ON_OFF;

  label = base_label + "/pulsatileFlow";
  if(!ptr->tp.getInspectedValue(label, ON_OFF)) {
    throw std::runtime_error(label + " is not set");
  }

  if(ON_OFF == "ON") {
    ptr->pulsatileFlow = ON;
  } else if(ON_OFF == "OFF") {
    ptr->pulsatileFlow = OFF;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  label = base_label + "/pulseBeginItr";
  if(!ptr->tp.getInspectedValue(label, ptr->pulseBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/T";
  if(!ptr->tp.getInspectedValue(label, ptr->T)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderInterface::readDarcyInfo()
{
  std::string str, base_label, label;
  base_label = "/DarcyParameter";

  label = base_label + "/alpha";
  if(!ptr->tp.getInspectedValue(label, ptr->alpha)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/resistance";
  if(!ptr->tp.getInspectedValue(label, ptr->resistance)) {
    throw std::runtime_error(label + " is not set");
  }
}

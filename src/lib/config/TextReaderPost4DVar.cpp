/**
 * @file TextReaderUSNS.cpp
 * @author K.Ueda
 * @date October, 2024
 */

#include "Config.h"

/**
 * @brief Read text parameter.
 */
void TextReaderPost4DVar::readBasicInfo(Config &conf)
{
  TextReaderInterface::readBasicInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReaderPost4DVar::readGridInfo(Config &conf)
{
  TextReaderInterface::readGridInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReaderPost4DVar::readResultsInfo(Config &conf)
{
  std::string str, base_label, sub_label, label;

  base_label = "/Results";

  label = base_label + "/inputDir";
  if(!conf.tp.getInspectedValue(label, conf.inputDir)) {
    throw std::runtime_error(label + " is not set");
  }
}

/**
 * @file   Tool.h
 * @author K.Ueda
 * @date   Jun, 2024
*/

#ifndef TOOL_H
#define TOOL_H

#include <iostream>
#include <vector>

class VecTool
{
public:
    template <typename T>
    static void resize(std::vector<T>& vec, int newSize){
        vec.resize(newSize, 0);
    }
    template <typename T>
    static void resize(std::vector<std::vector<T>>& vec, 
                       int newSizeOuter, int newSizeInner)
    {
        vec.resize(newSizeOuter);
        for(auto& innerVec : vec) {
            innerVec.resize(newSizeInner, 0);
        }
    }
    template <typename T>
    static void resize(std::vector<std::vector<std::vector<T>>>& vec, 
                       int newSizeOuter, int newSizeMiddle, int newSizeInner)
    {
        vec.resize(newSizeOuter);
        for(auto& middleVec : vec){
            middleVec.resize(newSizeMiddle);
            for(auto& innerVec : middleVec){
                innerVec.resize(newSizeInner, 0);
            }
        }
    }

    template <typename T>
    static void resize(std::vector<std::vector<std::vector<std::vector<T>>>>& vec, 
                         int newSizeOuter, int newSizeMiddle, int newSizeInner, int newSizeInnermost)
    {
        vec.resize(newSizeOuter);
        for(auto& middleVec : vec){
            middleVec.resize(newSizeMiddle);
            for(auto& innerVec : middleVec){
                innerVec.resize(newSizeInner);
                for(auto& innermostVec : innerVec){
                    innermostVec.resize(newSizeInnermost);
                }
            }
        }
    }
};

#endif

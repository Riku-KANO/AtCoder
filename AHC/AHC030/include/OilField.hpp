#pragma once

/**
 * @file OilField.hpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */

#include <global.hpp>
#include <geometry.hpp>

struct OilField
{
    Vec<Point> positions;
    int width;
    int height;

    OilField(): width(0), height(0){
        
    }
};
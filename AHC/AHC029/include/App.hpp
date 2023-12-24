#pragma once

#include "global.hpp"
#include "Property.hpp"
#include "input.hpp"
#include "Game.hpp"

class App
{
public:
    App(const Property& _props);
    void init();
    void run();
    void summary();

private:
    const Property props;
    Game game;
};
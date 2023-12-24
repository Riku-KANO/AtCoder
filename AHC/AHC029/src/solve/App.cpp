#include "App.hpp"

App::App(const Property& _props): props(_props){}

void App::init() {
    
    Input input = Input::get_input();
    
    game.init(input, props);

    LOG_INFO("Application Initialization done at %6.4f[s]", get_time(start_time));
}

void App::run() {
    
    game.run();

    LOG_INFO("Application Run Successfuly Ended at %6.4f[s]", get_time(start_time));
}

void App::summary() {

}
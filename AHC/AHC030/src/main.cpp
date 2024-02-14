#include <App.hpp>
#include <Property.hpp>

/**
 * # Idea
 * maximum likelihood estimation
 * uncertainty quantification
 * gradient descent (hill climbing)
 * 
 * # goal
 * - minimize entropy
 * - 
 */
int main(int argc, char* argv[]) {
    const Property props = Property::read_property(argc, argv);
    const Input input = Input::read_input();

    App app(input, props);
    app.init();
    app.run();
    return 0;
}

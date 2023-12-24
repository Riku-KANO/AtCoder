#include "App.hpp"

// ******************************* main ************************************
int main(int argc, char *argv[])
{
    Property props;
    start_time = clock();

#ifdef LOCAL
    props = arg_parse(argc, argv);
#endif

    App app(props);

    app.init();
    app.run();
#ifdef LOCAL
    app.summary();
#endif
    return 0;
}
// *************************************************************************
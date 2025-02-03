#include "choc_WebView.h"
#include "choc_DesktopWindow.h"
#include "choc_MessageLoop.h"
#include <iostream>

int main() {
    // Create the window first with specific bounds
    choc::ui::Bounds windowBounds{ 100, 100, 800, 600 }; // x, y, width, height
    choc::ui::DesktopWindow window(windowBounds);

    // Create webview with debug enabled
    choc::ui::WebView::Options options;
    options.enableDebugMode = true;
    choc::ui::WebView webview(options);

    if (!webview.loadedOK()) {
        std::cerr << "Failed to initialize WebView" << std::endl;
        return 1;
    }

    // Set webview as the window's content
    window.setContent(webview.getViewHandle());

    std::string html =
        "<html><body style='background:black; color:white'>"
        "<h1>DAWGE Test</h1>"
        "<p>If you can see this, the WebView is working!</p>"
        "</body></html>";

    webview.setHTML(html);

    // Run the message loop
    choc::messageloop::run();

    return 0;
}

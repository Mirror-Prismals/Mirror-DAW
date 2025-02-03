#define _CRT_SECURE_NO_WARNINGS
#include "choc_WebView.h"
#include "choc_DesktopWindow.h"
#include "choc_MessageLoop.h"
#include "choc_JSON.h"  // For JSON functionality
#include <iostream>

int main() {
    choc::ui::Bounds windowBounds{ 100, 100, 800, 600 };
    choc::ui::DesktopWindow window(windowBounds);

    choc::ui::WebView::Options options;
    options.enableDebugMode = true;

    choc::ui::WebView webview(options);
    if (!webview.loadedOK()) {
        std::cerr << "Failed to initialize WebView" << std::endl;
        return 1;
    }

    std::string html = R"(
<html><head><base href="https://mirrorguid.ia/generate/VGhlIGd1aS">
    <title>MirrorOS 'FF' - Emerald Cut Ruby Gemstone Interface</title>
    <style>
        body, html {
            margin: 0;
            padding: 0;
            width: 100%;
            height: 100%;
            overflow: hidden;
            background: #e0e0e0;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            transition: background-color 2s ease;
        }
        .ruby-interface {
            width: 100vw;
            height: 100vh;
            background: #e0e0e0;
            position: relative;
            overflow: hidden;
            display: flex;
            justify-content: center;
            align-items: center;
            perspective: 1000px;
        }
        .ruby-gem {
            width: 100vw;
            height: calc(100vw * 9 / 21);
            max-height: 100vh;
            background: linear-gradient(145deg, 
                rgba(200, 0, 0, 0.9), 
                rgba(150, 0, 0, 0.95),
                rgba(100, 0, 0, 1)
            );
            position: relative;
            clip-path: polygon(
                10% 0%, 90% 0%, 100% 10%, 100% 90%, 90% 100%, 10% 100%, 0% 90%, 0% 10%
            );
            box-shadow: 
                0 0 50px rgba(255, 0, 0, 0.3),
                inset 0 0 30px rgba(255, 255, 255, 0.1),
                inset 0 0 10px rgba(0, 0, 0, 0.3);
            display: flex;
            justify-content: center;
            align-items: center;
            transform-style: preserve-3d;
            transition: transform 2s cubic-bezier(0.25, 0.1, 0.25, 1);
        }
        .tunnel {
            width: 60%;
            height: 60%;
            background: linear-gradient(to center, rgba(200, 0, 0, 0.7), rgba(240, 240, 240, 0.9));
            box-shadow: inset 0 0 20px rgba(0, 0, 0, 0.2);
            position: relative;
            overflow: hidden;
        }
        .level {
            position: absolute;
            border: 2px solid;
            box-shadow: inset 0 0 10px rgba(200, 0, 0, 0.3);
            transition: all 0.3s ease;
        }
        .level-1 { 
            top: 5%; left: 5%; right: 5%; bottom: 5%; 
            border-color: rgba(200, 0, 0, 0.9);
            background: linear-gradient(to center, rgba(200, 0, 0, 0.6), rgba(200, 0, 0, 0.1));
        }
        .level-2 { 
            top: 10%; left: 10%; right: 10%; bottom: 10%; 
            border-color: rgba(200, 0, 0, 0.7);
            background: linear-gradient(to center, rgba(200, 0, 0, 0.5), rgba(200, 0, 0, 0.05));
        }
        .level-3 { 
            top: 15%; left: 15%; right: 15%; bottom: 15%; 
            border-color: rgba(200, 0, 0, 0.5);
            background: linear-gradient(to center, rgba(200, 0, 0, 0.4), rgba(200, 0, 0, 0.02));
        }
        .level-4 { 
            top: 20%; left: 20%; right: 20%; bottom: 20%; 
            border-color: rgba(200, 0, 0, 0.3);
            background: linear-gradient(to center, rgba(200, 0, 0, 0.3), rgba(200, 0, 0, 0.01));
        }
        .level-5 { 
            top: 25%; left: 25%; right: 25%; bottom: 25%; 
            border-color: rgba(200, 0, 0, 0.1);
            background: #e0e0e0;
            border: none;
            box-shadow: none;
            cursor: pointer;
        }
        .facet {
            position: absolute;
            background: linear-gradient(145deg, 
                rgba(255, 255, 255, 0.1), 
                rgba(255, 255, 255, 0.05),
                rgba(255, 255, 255, 0.02)
            );
        }
        .facet-1 { top: 0; left: 10%; right: 10%; height: 10%; transform: skew(-30deg); }
        .facet-2 { bottom: 0; left: 10%; right: 10%; height: 10%; transform: skew(30deg); }
        .facet-3 { top: 10%; bottom: 10%; left: 0; width: 10%; transform: skew(0, 30deg); }
        .facet-4 { top: 10%; bottom: 10%; right: 0; width: 10%; transform: skew(0, -30deg); }
        .gleam {
            position: absolute;
            width: 150%;
            height: 150%;
            background: radial-gradient(
                ellipse at center,
                rgba(255, 255, 255, 0.2) 0%,
                rgba(255, 255, 255, 0.1) 25%,
                rgba(255, 255, 255, 0.05) 50%,
                transparent 70%
            );
            transform: translateX(-50%) translateY(-50%);
            opacity: 0;
            transition: opacity 0.5s ease;
            pointer-events: none;
        }
        .info-panel {
            position: absolute;
            bottom: 20px;
            left: 20px;
            background: rgba(0, 0, 0, 0.7);
            color: white;
            padding: 15px;
            border-radius: 10px;
            font-size: 1.2vw;
            opacity: 0;
            transition: opacity 0.5s ease;
        }
        #continue-button, #zoom-out-button {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            padding: 10px 20px;
            font-size: 18px;
            background-color: #ff4444;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            opacity: 0;
            transition: opacity 0.5s ease;
            display: none;
        }
    </style>
</head>
<body>
    <div class="ruby-interface">
        <div class="ruby-gem">
            <div class="facet facet-1"></div>
            <div class="facet facet-2"></div>
            <div class="facet facet-3"></div>
            <div class="facet facet-4"></div>
            <div class="tunnel">
                <div class="level level-1"></div>
                <div class="level level-2"></div>
                <div class="level level-3"></div>
                <div class="level level-4"></div>
                <div class="level level-5"></div>
            </div>
            <div class="gleam"></div>
        </div>
        <div class="info-panel">
            Gemstone: Ruby<br>
            Cut: Emerald Cut<br>
            Color: Deep Red<br>
            Clarity: VVS1<br>
            Carat: 5.23<br>
            Origin: Myanmar (Burma)
        </div>
        <button id="continue-button">Continue</button>
        <button id="zoom-out-button">Zoom Out</button>
    </div>

    <script>
        const gem = document.querySelector('.ruby-gem');
        const gleam = document.querySelector('.gleam');
        const infoPanel = document.querySelector('.info-panel');
        const levels = document.querySelectorAll('.level');
        const centralLevel = document.querySelector('.level-5');
        const continueButton = document.getElementById('continue-button');
        const zoomOutButton = document.getElementById('zoom-out-button');
        let isZoomedIn = false;

        function updateGleam(event) {
            const rect = gem.getBoundingClientRect();
            const x = event.clientX - rect.left;
            const y = event.clientY - rect.top;
            
            gleam.style.left = x + 'px';
            gleam.style.top = y + 'px';
            gleam.style.opacity = '1';
            
            setTimeout(() => {
                gleam.style.opacity = '0';
            }, 500);
        }

        function updateLevels(event) {
            const rect = gem.getBoundingClientRect();
            const x = (event.clientX - rect.left) / rect.width;
            const y = (event.clientY - rect.top) / rect.height;

            levels.forEach((level, index) => {
                const scale = 1 + (index * 0.05 * (x + y - 1));
                level.style.transform = `scale(${scale})`;
                level.style.opacity = 1 - (index * 0.15);
            });
        }

        gem.addEventListener('mousemove', (event) => {
            if (!isZoomedIn) {
                updateGleam(event);
                updateLevels(event);
            }
        });

        gem.addEventListener('mouseleave', () => {
            if (!isZoomedIn) {
                gleam.style.opacity = '0';
                levels.forEach((level, index) => {
                    level.style.transform = 'scale(1)';
                    level.style.opacity = 1 - (index * 0.15);
                });
            }
        });

        function zoomIn() {
            isZoomedIn = true;
            gem.style.transform = 'translateZ(1000px) rotateX(30deg)';
            
            levels.forEach((level, index) => {
                if (index < 4) {
                    level.style.opacity = '0';
                    level.style.transition = 'opacity 1.5s ease';
                }
            });

            centralLevel.style.transition = 'all 2s cubic-bezier(0.25, 0.1, 0.25, 1)';
            centralLevel.style.transform = 'scale(5)';
            centralLevel.style.opacity = '1';

            setTimeout(() => {
                centralLevel.style.transform = 'scale(20)';
                centralLevel.style.opacity = '0';
            }, 1000);

            setTimeout(() => {
                document.body.style.backgroundColor = '#ffffff';
                continueButton.style.display = 'block';
                continueButton.style.opacity = '1';
                zoomOutButton.style.display = 'block';
                zoomOutButton.style.opacity = '1';
            }, 2000);
        }

        function zoomOut() {
            isZoomedIn = false;
            gem.style.transform = 'translateZ(0) rotateX(0)';
            document.body.style.backgroundColor = '#e0e0e0';

            levels.forEach((level, index) => {
                level.style.opacity = 1 - (index * 0.15);
                level.style.transition = 'opacity 1.5s ease';
            });

            centralLevel.style.transition = 'all 2s cubic-bezier(0.25, 0.1, 0.25, 1)';
            centralLevel.style.transform = 'scale(1)';
            centralLevel.style.opacity = '1';

            continueButton.style.display = 'none';
            continueButton.style.opacity = '0';
            zoomOutButton.style.display = 'none';
            zoomOutButton.style.opacity = '0';
        }

        centralLevel.addEventListener('click', () => {
            if (!isZoomedIn) {
                zoomIn();
            }
        });

        gem.addEventListener('click', (event) => {
            if (event.target !== centralLevel) {
                infoPanel.style.opacity = infoPanel.style.opacity === '1' ? '0' : '1';
            }
        });

        zoomOutButton.addEventListener('click', zoomOut);

        continueButton.addEventListener('click', () => {
            // Add functionality for the continue button here
            console.log('Continue button clicked');
        });
    </script>
</body></html>
    )";

    webview.setHTML(html);
    window.setContent(webview.getViewHandle());

    // Create test data using Choc's JSON functionality
    auto testData = choc::json::parse(R"({
        "message": "Hello from C++",
        "number": 42
    })");

    // Send data to JavaScript
    webview.evaluateJavascript(
        "window.updateOutput?.(" + choc::json::toString(testData, true) + ");",
        [](const std::string& error, const choc::value::ValueView& result) {
            if (error.empty())
                std::cout << "JavaScript executed: " << result.toString() << std::endl;
            else
                std::cerr << "Error: " << error << std::endl;
        }
    );

    choc::messageloop::run();
    return 0;
}

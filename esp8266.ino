#include <Grandeur.h>
#include "WiFi.h"

// WiFi credentials
const char* ssid = "FRITZ!Box 7530 UQ";
const char* passphrase = "10860176429693588533";

// Grandeur credentials
const char * apiKey = "grandeurl47029t311w70pxf0eh12gie";
const char* token = "eyJ0b2tlbiI6ImV5SmhiR2NpT2lKSVV6STFOaUlzSW5SNWNDSTZJa3BYVkNKOS5leUpwWkNJNkltUmxkbWxqWld3ME56QjBhbko1TVRGNE5qQndlR1l6WTNWelpUVndZaUlzSW5SNWNHVWlPaUprWlhacFkyVWlMQ0pwWVhRaU9qRTJOVFEzTnpreU1UZDkuZVFsejRIYzl1aUstTE1uWVRMclVMVlZPTlF5bFFiV094THlTMTVSbFVIWSJ9";
const char* deviceId = "devicel470tjpx11x50pxf0mndhitx";

Grandeur::Project project;

void setup() {
    Serial.begin(9600);
    // This connects the device to WiFi.
    connectToWiFi(ssid, passphrase);
    // This connects the device to internet.
    project = grandeur.init(apiKey, token);
}

unsigned long current = millis();

void loop() {
    // This sends data to internet.
    if(project.isConnected() && millis() - current > 500) {
        project.device(deviceId).data().set("millis", millis());
        current = millis();
    }

    // This runs the SDK when the device WiFi is connected.
    if(WiFi.status() == WL_CONNECTED) project.loop();
}

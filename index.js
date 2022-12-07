const apiKey = "grandeurl47029t311w70pxf0eh12gie";
const accessKey = "eyJ0b2tlbiI6ImV5SmhiR2NpT2lKSVV6STFOaUlzSW5SNWNDSTZJa3BYVkNKOS5leUpwWkNJNkltRmpZMlZ6YzJ3ME56Rm5iM1I0TVRGNGRqQndlR1prYjNnNVltSXplaUlzSW5SNWNHVWlPaUpoWTJObGMzTWlMQ0pwWVhRaU9qRTJOVFEzT0RBeU9UWjkudTZ2YUh4c2hQUk1LQlg0eC1QdDF4Z3lNZWhsZUtjMXJtaTZJWHdlTkw1ayJ9";
const accessToken = "eyJ0b2tlbiI6ImV5SmhiR2NpT2lKSVV6STFOaUlzSW5SNWNDSTZJa3BYVkNKOS5leUpwWkNJNkltUmxkbWxqWld3ME56QjBhbko1TVRGNE5qQndlR1l6WTNWelpUVndZaUlzSW5SNWNHVWlPaUprWlhacFkyVWlMQ0pwWVhRaU9qRTJOVFEzTnpreU1UZDkuZVFsejRIYzl1aUstTE1uWVRMclVMVlZPTlF5bFFiV094THlTMTVSbFVIWSJ9";

// This connects the webpage to the internet.
const project = grandeur.init(apiKey, accessKey, accessToken);
project.auth().login("lhs75@iclouc.com", "Leopold75");

// This subscribes to the "millis" variable.
project.devices().device(devicel470tjpx11x50pxf0mndhitx).data().on("millis", (path, value) => document.getElementById("data").innerHTML = value);

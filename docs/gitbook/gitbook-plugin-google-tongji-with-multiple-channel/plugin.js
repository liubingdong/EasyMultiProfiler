require(["gitbook"], function(gitbook) {
    gitbook.events.bind("start", function(e, config) {
        // Get plugin config
        var config = config["google-tongji-with-multiple-channel"] || {};
        var token = config.token;
        var url = config.url;
        var multipleChannelConfig = config.multipleChannelConfig;

        // Enable multi-channel configuration
        if (JSON.stringify(multipleChannelConfig) !== "{}") {
            var currentHost = location.hostname;
            for(var channelHost in multipleChannelConfig){
                if (currentHost == channelHost) {
                    token = multipleChannelConfig[channelHost].token || token;
                    url = multipleChannelConfig[channelHost].url || url;
                    break;
                }
            }
        }
        
        // Automatically insert Google analysis code
        if(url !== "" && token !== ""){
            // Reference js
            var hm = document.createElement("script");
            hm.async = true;
            hm.src = url + '?id=' + token;
            var s = document.getElementsByTagName("script")[0];
            s.parentNode.insertBefore(hm, s);

            // Use js
            var usejs = document.createElement("script");
            var useData = [];
            useData.push("window.dataLayer = window.dataLayer || [];");
            useData.push("function gtag(){dataLayer.push(arguments);}");
            useData.push("gtag('js', new Date());");
            useData.push("gtag('config', '"+ token + "');");
            usejs.innerHTML = useData.join("\n");
            var gtag = document.getElementsByTagName("script")[1];
            gtag.parentNode.insertBefore(usejs, gtag);
        }
    });
});
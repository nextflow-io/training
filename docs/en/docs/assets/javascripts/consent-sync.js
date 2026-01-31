// Sync MkDocs Material settings across all language paths when written
// Intercepts localStorage.setItem to copy settings (consent, palette) to all languages
(function () {
  // Get language paths from the language picker in the DOM
  var langPaths = ["/"];
  document
    .querySelectorAll(".md-select__link[hreflang]")
    .forEach(function (link) {
      var lang = link.getAttribute("hreflang");
      if (lang && lang !== "en") {
        langPaths.push("/" + lang + "/");
      }
    });

  // Store original setItem
  var originalSetItem = localStorage.setItem.bind(localStorage);

  // Override setItem to sync MkDocs Material settings
  localStorage.setItem = function (key, value) {
    // Call original first
    originalSetItem(key, value);

    // Check if this is a MkDocs Material setting (e.g., /.__consent, /de/.__palette)
    var match = key.match(/^(\/[a-z]{2}\/|\/)(.__\w+)$/);
    if (match) {
      var suffix = match[2]; // e.g., ".__consent"
      // Copy to all other language paths
      langPaths.forEach(function (path) {
        var otherKey = path + suffix;
        if (otherKey !== key) {
          originalSetItem(otherKey, value);
        }
      });
    }
  };
})();

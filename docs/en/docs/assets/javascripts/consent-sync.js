// Sync cookie consent and palette (light/dark mode) across language paths
// MkDocs Material prefixes localStorage keys with path (e.g., /de/.__consent, /de/.__palette)
// This script copies settings from any language to the current one on page load
(function () {
  var suffixes = [".__consent", ".__palette"];

  // Get current path's prefix
  var pathMatch = window.location.pathname.match(/^\/[a-z]{2}\//);
  var currentPath = pathMatch ? pathMatch[0] : "/";

  suffixes.forEach(function (suffix) {
    var currentKey = currentPath + suffix;

    // Skip if current path already has this setting
    if (localStorage.getItem(currentKey)) {
      return;
    }

    // Find any existing setting in localStorage
    for (var i = 0; i < localStorage.length; i++) {
      var key = localStorage.key(i);
      if (key && key.endsWith(suffix)) {
        var value = localStorage.getItem(key);
        if (value) {
          localStorage.setItem(currentKey, value);
          break;
        }
      }
    }
  });
})();

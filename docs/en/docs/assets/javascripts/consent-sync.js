// Sync MkDocs Material settings (consent, palette, etc.) across language paths
// MkDocs Material prefixes localStorage keys with path (e.g., /de/.__consent, /de/.__palette)
// This script copies all settings from any language to the current one on page load
(function () {
  // Get current path's prefix
  var pathMatch = window.location.pathname.match(/^\/[a-z]{2}\//);
  var currentPath = pathMatch ? pathMatch[0] : "/";

  // Find all MkDocs Material settings (keys ending with .__something)
  var settingsToSync = {};
  for (var i = 0; i < localStorage.length; i++) {
    var key = localStorage.key(i);
    if (key) {
      var match = key.match(/\.__\w+$/);
      if (match) {
        var suffix = match[0];
        // Only store first found value for each suffix
        if (!settingsToSync[suffix]) {
          settingsToSync[suffix] = localStorage.getItem(key);
        }
      }
    }
  }

  // Copy any missing settings to current path
  Object.keys(settingsToSync).forEach(function (suffix) {
    var currentKey = currentPath + suffix;
    if (!localStorage.getItem(currentKey)) {
      localStorage.setItem(currentKey, settingsToSync[suffix]);
    }
  });
})();

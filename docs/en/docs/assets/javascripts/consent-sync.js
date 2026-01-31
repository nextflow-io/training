// Sync cookie consent across language paths
// MkDocs Material prefixes localStorage keys with path (e.g., /de/.__consent)
// This script copies consent from any language to all others
(function () {
  // Find any existing consent in localStorage
  for (var i = 0; i < localStorage.length; i++) {
    var key = localStorage.key(i);
    if (key && key.endsWith(".__consent")) {
      var consent = localStorage.getItem(key);
      // Get current path's consent key
      var pathMatch = window.location.pathname.match(/^\/[a-z]{2}\//);
      var currentPath = pathMatch ? pathMatch[0] : "/";
      var currentKey = currentPath + ".__consent";
      // If current path doesn't have consent, copy it
      if (!localStorage.getItem(currentKey) && consent) {
        localStorage.setItem(currentKey, consent);
      }
      break;
    }
  }
})();

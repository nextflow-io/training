// Fix language picker links for mike versioning
// Detects version prefix (e.g., /latest/, /0.dev/) and prepends to language links
(function () {
  function fixLanguageLinks() {
    var path = window.location.pathname;
    var versionMatch = path.match(/^\/([^\/]+)\//);
    var knownLangs = ["en", "pt", "es", "fr", "it", "ko", "pl", "tr"];
    var versionPrefix = "";

    if (versionMatch) {
      var potentialVersion = versionMatch[1];
      // If first path segment isn't a language code, it's a version
      if (knownLangs.indexOf(potentialVersion) === -1) {
        versionPrefix = "/" + potentialVersion;
      }
    }

    // Only proceed if we have a version prefix to add
    if (!versionPrefix) return;

    var langLinks = document.querySelectorAll(".md-select__link[hreflang]");
    langLinks.forEach(function (link) {
      var href = link.getAttribute("href");
      // Handle both absolute URLs and relative paths
      try {
        var url = new URL(href, window.location.origin);
        var linkPath = url.pathname;

        // Check if this link already has the version prefix
        if (
          !linkPath.startsWith(versionPrefix + "/") &&
          linkPath !== versionPrefix
        ) {
          // Prepend version prefix to the path
          url.pathname = versionPrefix + linkPath;
          link.setAttribute("href", url.href);
        }
      } catch (e) {
        // Fallback for older browsers
        if (href.startsWith("/") && !href.startsWith(versionPrefix)) {
          link.setAttribute("href", versionPrefix + href);
        }
      }
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", fixLanguageLinks);
  } else {
    fixLanguageLinks();
  }

  // Re-run after instant navigation (mkdocs-material uses RxJS)
  if (typeof document$ !== "undefined") {
    document$.subscribe(fixLanguageLinks);
  }
})();

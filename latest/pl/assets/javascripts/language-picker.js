// Fix language picker links to preserve current page when switching languages
// Also handles mike versioning prefixes (e.g., /latest/, /0.dev/)
(function () {
  function fixLanguageLinks() {
    var path = window.location.pathname;

    // Get known languages from the language picker in the DOM
    var knownLangs = [];
    var langLinks = document.querySelectorAll(".md-select__link[hreflang]");
    langLinks.forEach(function (link) {
      var lang = link.getAttribute("hreflang");
      if (lang && knownLangs.indexOf(lang) === -1) {
        knownLangs.push(lang);
      }
    });

    if (knownLangs.length === 0) {
      console.error("language-picker.js: Could not find language links in DOM");
      return;
    }

    // Check if a segment looks like a mike version
    // Versions are like: latest, stable, 0.dev, 1.0, 2.8.1, etc.
    function isVersion(segment) {
      return /^(latest|stable|\d+(\.\d+)*\.?dev|\d+(\.\d+)*)$/.test(segment);
    }

    // Parse the current path to extract version prefix, language, and page path
    // Handles: /page, /de/page, /latest/page, /latest/de/page
    var versionPrefix = "";
    var pagePath = path;

    // Split path into segments
    var segments = path.split("/").filter(function (s) {
      return s.length > 0;
    });

    if (segments.length > 0) {
      // Check if first segment is a version (not a language, and looks like a version)
      if (knownLangs.indexOf(segments[0]) === -1 && isVersion(segments[0])) {
        // First segment is a version prefix
        versionPrefix = "/" + segments[0];
        segments.shift();
      }

      // Check if first/next segment is a language
      if (segments.length > 0 && knownLangs.indexOf(segments[0]) !== -1) {
        // Remove language from page path
        segments.shift();
      }

      // Remaining segments are the page path
      pagePath = "/" + segments.join("/");
      // Preserve trailing slash if original had one
      if (path.endsWith("/") && !pagePath.endsWith("/")) {
        pagePath += "/";
      }
    }

    // Update each language link to preserve current page
    langLinks.forEach(function (link) {
      var lang = link.getAttribute("hreflang");
      var langPrefix = lang === "en" ? "" : "/" + lang;
      var newHref = versionPrefix + langPrefix + pagePath;
      // Ensure we don't end up with empty href
      if (newHref === "") {
        newHref = "/";
      }
      link.setAttribute("href", newHref);
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

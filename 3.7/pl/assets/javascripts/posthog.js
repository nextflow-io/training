// PostHog Analytics Integration

(function () {
  !(function (t, e) {
    var o, n, p, r;
    e.__SV ||
      ((window.posthog = e),
      (e._i = []),
      (e.init = function (i, s, a) {
        function g(t, e) {
          var o = e.split(".");
          2 == o.length && ((t = t[o[0]]), (e = o[1])),
            (t[e] = function () {
              t.push([e].concat(Array.prototype.slice.call(arguments, 0)));
            });
        }
        ((p = t.createElement("script")).type = "text/javascript"),
          (p.async = !0),
          (p.src =
            s.api_host.replace(".i.posthog.com", "-assets.i.posthog.com") +
            "/static/array.js"),
          (r = t.getElementsByTagName("script")[0]).parentNode.insertBefore(
            p,
            r,
          );
        var u = e;
        for (
          void 0 !== a ? (u = e[a] = []) : (a = "posthog"),
            u.people = u.people || [],
            u.toString = function (t) {
              var e = "posthog";
              return (
                "posthog" !== a && (e += "." + a), t || (e += " (stub)"), e
              );
            },
            u.people.toString = function () {
              return u.toString(1) + ".people (stub)";
            },
            o =
              "capture identify alias people.set people.set_once set_config register register_once unregister opt_out_capturing has_opted_out_capturing opt_in_capturing reset isFeatureEnabled onFeatureFlags getFeatureFlag getFeatureFlagPayload reloadFeatureFlags group updateEarlyAccessFeatureEnrollment getEarlyAccessFeatures getActiveMatchingSurveys getSurveys getNextSurveyStep onSessionId".split(
                " ",
              ),
            n = 0;
          n < o.length;
          n++
        )
          g(u, o[n]);
        e._i.push([i, s, a]);
      }),
      (e.__SV = 1));
  })(document, window.posthog || []);

  var isPostHogInitialized = false;

  function initializePostHog() {
    if (isPostHogInitialized) return;

    // Check if Material for MkDocs consent cookie exists
    var consent = __md_get && __md_get("__consent");

    if (consent && consent.posthog) {
      // User has consented to PostHog specifically - use your exact configuration
      posthog.init("phc_pLVMvF3StGtifZ71oQrnQ99FGD1luzjIcG3xrs1PIaW", {
        api_host: "https://eu.i.posthog.com",
        person_profiles: "always",
        persistence: "localStorage",
        opt_out_capturing_by_default: false,
        disable_web_experiments: false,
      });
    } else {
      // User hasn't consented to PostHog or declined - use memory mode
      posthog.init("phc_pLVMvF3StGtifZ71oQrnQ99FGD1luzjIcG3xrs1PIaW", {
        api_host: "https://eu.i.posthog.com",
        person_profiles: "always",
        persistence: "memory",
        opt_out_capturing_by_default: true,
        disable_web_experiments: false,
      });
    }

    isPostHogInitialized = true;
  }

  // Function to handle consent changes after initialization
  function handlePostHogConsent() {
    if (!isPostHogInitialized) {
      initializePostHog();
      return;
    }

    // Check PostHog-specific consent and update config
    var consent = __md_get && __md_get("__consent");

    if (consent && consent.posthog) {
      // User consented to PostHog specifically - switch to full tracking
      posthog.set_config({
        persistence: "localStorage",
      });
      posthog.opt_in_capturing();
    } else {
      // User declined PostHog or no consent - switch to memory mode
      posthog.set_config({
        persistence: "memory",
      });
      posthog.opt_out_capturing();
    }
  }

  // Wait for MkDocs to load, then initialize and handle consent
  document.addEventListener("DOMContentLoaded", function () {
    // Initial initialization
    initializePostHog();

    // Listen for consent changes (when user interacts with cookie banner)
    document.addEventListener("click", function (e) {
      // Check if it's a consent button or PostHog checkbox
      if (
        e.target.hasAttribute("data-md-consent") ||
        (e.target.type === "checkbox" && e.target.name === "posthog")
      ) {
        setTimeout(handlePostHogConsent, 100);
      }
    });

    // Also listen for any changes in consent (alternative method)
    document.addEventListener("change", function (e) {
      // Check if it's the PostHog checkbox specifically
      if (e.target.type === "checkbox" && e.target.name === "posthog") {
        setTimeout(handlePostHogConsent, 100);
      }
    });
  });
})();

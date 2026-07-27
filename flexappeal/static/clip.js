/* Loads the decorative trajectory clip, last and only if it is wanted.
 *
 * The <video> in _clip.html ships with a poster and no src, so the browser
 * renders the still immediately and fetches nothing. This attaches the real
 * source only after window.load, which means the 811 KB cannot compete with
 * the page's own CSS, JavaScript or (on the Analysis tab) the several-megabyte
 * Mol* and Plotly bundles for bandwidth during first paint.
 *
 * Two ways it deliberately does nothing at all: a visitor who has asked for
 * reduced motion never downloads the video, and with JavaScript off the poster
 * simply stays put. Both leave the same still image the video opens on, so
 * nothing is missing in either case -- the panel repeats no information that
 * appears nowhere else.
 */
(function () {
  'use strict';

  const still = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  function load() {
    document.querySelectorAll('video.md-clip[data-src]').forEach(function (video) {
      video.src = video.dataset.src;
      video.removeAttribute('data-src');
      // Safari returns a rejected promise rather than throwing when autoplay is
      // refused; without the catch that surfaces as an unhandled rejection in
      // the console on a page where nothing has actually gone wrong.
      const started = video.play();
      if (started && started.catch) started.catch(function () {});
    });
  }

  if (still) return;
  if (document.readyState === 'complete') load();
  else window.addEventListener('load', load);
})();

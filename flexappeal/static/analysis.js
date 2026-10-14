/* Analysis tab: draw the Plotly panels, drive the Mol* viewer, and link them.
 *
 * Vanilla JS, no framework. Both libraries are vendored, never fetched from a
 * CDN, so this page works with no third-party network access at all.
 */

(function () {
  'use strict';

  /* Playback is a frame rate with a clamp on the pass, not a fixed pass.
     
     It used to be a flat fifty seconds however many frames the file held, which
     divides out to the interval between frames: 100 ms for the 500-frame
     example, and 1,000 ms for the 50-frame bundles the reference-state viewer
     carried. One frame per second is not slow playback, it is a slideshow, and
     it was reported as jerky and in slow motion because that is what it was.
     
     So the rate is the constant and the pass length follows from it. Twenty-five
     a second is smooth to the eye and well inside what a model-index update can
     sustain. The clamp is on the pass rather than on the rate: a short file
     played at 25 fps would be over in two seconds, and a very long one would
     take minutes, so the pass is held between eight and forty seconds and the
     rate gives way at the ends. Mol*'s own durationInS clamps to 1-120, which
     this is well inside. */
  const PLAYBACK_FPS = 25;
  const PLAYBACK_MIN_SECONDS = 8;
  const PLAYBACK_MAX_SECONDS = 40;

  /* How the rest of this file asks playback to stop and start again.
     
     Loading a structure into Mol* is main-thread work, and so is stepping the
     model index twenty-five times a second: with the trajectory running, a
     reference ghost took more than ten seconds to appear and its toggle sat
     disabled the whole time, which reads as a broken control. It was not
     noticed before because playback used to run at one frame a second, where
     there was plenty of thread to spare.
     
     Set by `wireSlider`, which owns the frame clock. Left as no-ops when there
     is no slider on the page, so a caller never has to check. */
  let playback = {pause: function () { return false; }, resume: function () {}};

  function playbackSeconds(frames) {
    if (!frames || frames < 2) return PLAYBACK_MIN_SECONDS;
    return Math.min(PLAYBACK_MAX_SECONDS,
                    Math.max(PLAYBACK_MIN_SECONDS, frames / PLAYBACK_FPS));
  }

  /* ------------------------------------------------------- drop zone */

  function initDropZone() {
    const zone = document.getElementById('drop-zone');
    const input = document.getElementById('fxa-input');
    const hint = document.getElementById('drop-hint');
    if (!zone || !input) return;

    input.addEventListener('change', function () {
      if (input.files.length) hint.textContent = input.files[0].name;
    });

    ['dragenter', 'dragover'].forEach(function (event) {
      zone.addEventListener(event, function (e) {
        e.preventDefault();
        zone.classList.add('is-over');
      });
    });
    ['dragleave', 'drop'].forEach(function (event) {
      zone.addEventListener(event, function (e) {
        e.preventDefault();
        zone.classList.remove('is-over');
      });
    });
    zone.addEventListener('drop', function (e) {
      if (!e.dataTransfer || !e.dataTransfer.files.length) return;
      input.files = e.dataTransfer.files;
      hint.textContent = e.dataTransfer.files[0].name;
    });
  }

  /* ---------------------------------------------------------- plots */

  const PLOT_CONFIG = {
    responsive: true,
    displaylogo: false,
    // lasso and box-select do nothing useful on a time series and only crowd
    // the mode bar; the rest stays for zoom, pan and PNG export.
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
    toImageButtonOptions: { format: 'png', scale: 2 }
  };

  const drawn = [];

  /* Drawn a few at a time, after the first paint, rather than all at once
     before it.
     
     `Plotly.newPlot` is synchronous and it is not cheap: a DSSP heatmap is
     tens of thousands of cells. Drawing every panel in one pass out of
     DOMContentLoaded holds the main thread until the last one is done, and the
     browser cannot paint while it is held. On the example's fourteen panels
     that cost about a second; on a reference-state run page, which carries
     thirty-five, first contentful paint measured 12 to 15 seconds against
     section 6.3's 2 s budget. The page is complete and readable before any
     chart exists, which is exactly the argument 4bp made for deferring the
     viewer.
     
     Two panels per frame, so the browser paints between batches and the page
     scrolls while the rest arrive. The order is the document's, so the panels
     a reader is looking at are the ones drawn first. */
  const PLOTS_PER_FRAME = 2;

  function drawPlot(node) {
    let figure;
    try {
      figure = JSON.parse(node.dataset.figure);
    } catch (err) {
      node.innerHTML = '<p class="md-plot-error">This panel could not be drawn.</p>';
      return;
    }
    /* A panel that asks for a height gets it.
       
       The stylesheet gives `.md-plot` a fixed height, which is right for a time
       series and wrong for anything whose height is its content: the
       protein-ligand contacts panel is a horizontal bar per residue and it
       already computed a height from the number of bars, which the container
       then overrode. Twenty bars were drawn into a box the size of ten. */
    if (figure.layout && figure.layout.height) {
      node.style.height = figure.layout.height + 'px';
    }
    Plotly.newPlot(node, figure.data, figure.layout, PLOT_CONFIG);
    drawn.push(node);
    wirePlotClick(node);
  }

  function initPlots() {
    if (typeof Plotly === 'undefined') return;

    const pending = Array.prototype.slice.call(
      document.querySelectorAll('.md-plot[data-figure]'));

    window.addEventListener('resize', function () {
      drawn.forEach(function (node) { Plotly.Plots.resize(node); });
    });

    if (typeof requestAnimationFrame !== 'function') {
      pending.forEach(drawPlot);
      return;
    }
    (function batch() {
      for (let i = 0; i < PLOTS_PER_FRAME && pending.length; i += 1) {
        drawPlot(pending.shift());
      }
      if (pending.length) requestAnimationFrame(batch);
    })();
  }

  /* Clicking a point that knows its frame moves the viewer to that frame.
     
     Section 6.2 asks for this on the reference state map, where the question is
     "what was the protein doing when it was over there": a reader picks a point
     off the cloud and wants to see that conformation. The panel that carries it
     says so by putting a frame index in each point's `customdata`, so this is
     wired to every panel and does nothing on the ones that do not.
     
     Driven through the slider rather than by reaching into the viewer, so one
     control owns the model index and the two cannot end up disagreeing about
     which frame is shown. The frame is clamped to what the slider actually has:
     a results file is packed at a stride, so a trajectory frame index is not
     necessarily a frame the viewer holds. */
  function wirePlotClick(node) {
    if (!node || typeof node.on !== 'function') return;
    node.on('plotly_click', function (event) {
      const point = (event.points || [])[0];
      if (!point || point.customdata === undefined
          || point.customdata === null) return;
      const frame = Array.isArray(point.customdata)
        ? Number(point.customdata[0]) : Number(point.customdata);
      if (!Number.isFinite(frame)) return;
      const slider = document.getElementById('scrub');
      if (!slider) return;
      const highest = parseInt(slider.max, 10);
      slider.value = String(Math.max(0, Math.min(highest, Math.round(frame))));
      slider.dispatchEvent(new Event('input', { bubbles: true }));
      const host = document.getElementById('molstar-viewer');
      if (host && typeof host.scrollIntoView === 'function') {
        host.scrollIntoView({ behavior: 'smooth', block: 'center' });
      }
    });
  }

  /* --------------------------------------------------------- viewer */

  /* Neon presentation, matching DeepFlex's pocket viewer so the same complex
   * looks like the same complex in both places.
   *
   * The constants are DeepFlex's, which got them from BINMAN, and every one of
   * them is tuning somebody already paid for:
   *
   *   - carbons are pinned to one vivid colour rather than coloured by chain. A
   *     chain-coloured ligand can come out the same hue as its own shell and
   *     vanish into the glow; the colour a ligand gets should not depend on
   *     which chain it happens to sit beside. Oxygen, nitrogen and sulfur keep
   *     their element colours, so the molecule still reads chemically.
   *   - sticks at 0.34 rather than the default. Thinner and the glow sits over
   *     them, and the ligand reads as a cloud of unconnected dots: the bonds are
   *     there and simply cannot be seen.
   *   - the halo is a shell, not a fog. The probe radius inflates the surface
   *     clear of the atoms, and its emissive is below the sticks' own so the
   *     bloom pass cannot wash them out.
   *   - bloom runs in `emissive` mode, so only geometry carrying emissive > 0
   *     blooms and the protein behind the ligand stays readable.
   *
   * All of it is guarded. A viewer that cannot apply the styling shows the
   * default presentation, which is a worse-looking panel rather than no panel.
   */
  const NEON = {
    background: 0x07080F,
    neon: 0x00F5D4,
    carbon: 0xFF2D9E,
    stick: { sizeFactor: 0.34, aspectRatio: 1, emissive: 0.6 },
    cloud: { alpha: 0.09, emissive: 0.7, probeRadius: 2.4, resolution: 0.5 },
    bloom: { strength: 1.3, radius: 0.6, threshold: 0,
             mode: 'emissive', transparency: true }
  };

  function ligandComponents(plugin) {
    const out = [];
    try {
      const current = plugin.managers.structure.hierarchy.current;
      (current.structures || []).forEach(function (ref) {
        (ref.components || []).forEach(function (component) {
          const key = String(component.key || '').toLowerCase();
          let label = '';
          try {
            label = String(component.cell.obj.label || '').toLowerCase();
          } catch (err) { /* unlabelled is still testable by key */ }
          if (key.indexOf('ligand') >= 0 || label.indexOf('ligand') >= 0) {
            out.push(component);
          }
        });
      });
    } catch (err) { /* no hierarchy yet */ }
    return out;
  }

  /* Section 6.2's cartoon coloured by motif.
   *
   * The number comes out of the structure, not out of an API call. The
   * vendored Mol* bundle is minified and exports no MolScript, so building a
   * component per motif is not reachable from a page; the colour theme
   * registry does expose `uncertainty`, which reads the B-factor column, and
   * `tools/analyse_run.py` writes each atom's `--motif-N` number there.
   *
   * A set palette over a domain of [-0.5, n+0.5] puts each integer in the
   * middle of its own bin, so the mapping from motif number to colour is exact
   * rather than an artefact of where a continuous scale happens to fall.
   *
   * No-ops on the Analysis tab, which has no motif numbers to colour by: the
   * bundle says whether it painted any, and an unpainted structure keeps the
   * default theme rather than being washed to a single colour.
   */
  const MOTIF_SLOTS = 6;

  /* Mol* wants 0xRRGGBB, not a CSS string. Handed "#0072b2" it produced a
     cartoon that was uniformly almost black, which is the worst outcome
     available: the theme applied, so nothing errored, and the structure was
     still there and unreadable. Anything this cannot parse returns null and
     the caller leaves the default theme alone. */
  function cssColourToInt(value) {
    const text = String(value || '').trim();
    let hex = null;
    if (/^#[0-9a-f]{6}$/i.test(text)) hex = text.slice(1);
    else if (/^#[0-9a-f]{3}$/i.test(text)) {
      hex = text.slice(1).split('').map(function (c) { return c + c; }).join('');
    } else {
      const rgb = text.match(/^rgba?\(\s*(\d+)[,\s]+(\d+)[,\s]+(\d+)/i);
      if (rgb) {
        return (Number(rgb[1]) << 16) + (Number(rgb[2]) << 8) + Number(rgb[3]);
      }
    }
    if (hex === null) return null;
    const parsed = parseInt(hex, 16);
    return Number.isFinite(parsed) ? parsed : null;
  }

  function motifPalette() {
    const root = getComputedStyle(document.documentElement);
    /* Index 0 is everything no motif covers: the muted grey the rest of the
       page uses for "nothing to say". */
    const first = cssColourToInt(root.getPropertyValue('--md-text-muted'));
    const colours = [first === null ? 0x5a6b80 : first];
    for (let i = 1; i <= MOTIF_SLOTS; i += 1) {
      const value = cssColourToInt(root.getPropertyValue('--motif-' + i));
      if (value === null) return null;
      colours.push(value);
    }
    return colours;
  }

  function polymerComponents(plugin) {
    const out = [];
    try {
      const current = plugin.managers.structure.hierarchy.current;
      (current.structures || []).forEach(function (ref) {
        (ref.components || []).forEach(function (component) {
          const key = String(component.key || '').toLowerCase();
          let label = '';
          try {
            label = String(component.cell.obj.label || '').toLowerCase();
          } catch (err) { /* unlabelled is still testable by key */ }
          if (key.indexOf('polymer') >= 0 || label.indexOf('polymer') >= 0) {
            out.push(component);
          }
        });
      });
    } catch (err) { /* no hierarchy yet */ }
    return out;
  }

  function colourByMotif(viewer, host) {
    const plugin = viewer && viewer.plugin;
    if (!plugin || !host || host.dataset.motifPainted !== '1') return false;
    const colours = motifPalette();
    if (!colours) return false;

    const components = polymerComponents(plugin);
    if (!components.length) return false;
    let touched = 0;
    try {
      const update = plugin.build();
      components.forEach(function (component) {
        (component.representations || []).forEach(function (rep) {
          if (!rep.cell) return;
          update.to(rep.cell).update(function (old) {
            if (!old.type || old.type.name !== 'cartoon') return;
            /* Reversed, because `uncertainty` puts high values at the start
               of the scale: a higher B-factor means less certainty, and the
               theme's palette runs from certain to uncertain. Unreversed, the
               background took motif 6's blue and each motif took the colour of
               its opposite number, which looked entirely deliberate. Confirmed
               by predicting the whole mapping from the atom counts and
               matching it against the render rather than by flipping it until
               it looked right.

               Compensated here rather than by writing `6 - n` into the PDB,
               so the column keeps the only meaning worth having: the number in
               the B-factor is the motif's number. */
            old.colorTheme = {
              name: 'uncertainty',
              params: {
                domain: [-0.5, MOTIF_SLOTS + 0.5],
                list: { kind: 'set', colors: colours.slice().reverse() }
              }
            };
            touched += 1;
          });
        });
      });
      update.commit();
    } catch (err) {
      return false;
    }
    return touched > 0;
  }

  function applyNeon(viewer) {
    const plugin = viewer.plugin;
    if (!plugin || !plugin.canvas3d) return;

    try {
      plugin.canvas3d.setProps({
        renderer: { backgroundColor: NEON.background },
        postprocessing: { bloom: { name: 'on', params: NEON.bloom } }
      });
    } catch (err) { /* a light viewport is not worth failing the mount over */ }

    const components = ligandComponents(plugin);
    if (!components.length) return;

    /* Retune the ligand's existing sticks in place. Only a ball-and-stick is
       touched: renaming a type while keeping another type's params is how a
       representation ends up with nonsense settings. */
    try {
      const update = plugin.build();
      components.forEach(function (component) {
        (component.representations || []).forEach(function (rep) {
          if (!rep.cell) return;
          update.to(rep.cell).update(function (old) {
            if (!old.type || old.type.name !== 'ball-and-stick') return;
            old.type.params = old.type.params || {};
            old.type.params.sizeFactor = NEON.stick.sizeFactor;
            old.type.params.aspectRatio = NEON.stick.aspectRatio;
            old.type.params.emissive = NEON.stick.emissive;
            /* No hydrogens on the ligand. They treble the stick count on a
               small molecule and bury the heavy-atom skeleton that the pose,
               the contacts and the 2D diagram beside it are all about. The
               diagram strips them too (measure.depiction uses RemoveHs), so
               the two views of the ligand now show the same atoms. */
            old.type.params.ignoreHydrogens = true;
            old.colorTheme = {
              name: 'element-symbol',
              params: {
                carbonColor: { name: 'uniform', params: { value: NEON.carbon } }
              }
            };
          });
        });
      });
      update.commit();
    } catch (err) { /* default sticks are a fine outcome */ }

    /* The shell. Added as a second representation rather than replacing the
       sticks, because it is the glow around them. */
    try {
      const builder = plugin.builders.structure.representation;
      components.forEach(function (component) {
        builder.addRepresentation(component.cell, {
          type: 'molecular-surface',
          typeParams: {
            alpha: NEON.cloud.alpha, emissive: NEON.cloud.emissive,
            probeRadius: NEON.cloud.probeRadius,
            resolution: NEON.cloud.resolution, ignoreLight: true
          },
          color: 'uniform',
          colorParams: { value: NEON.neon }
        }).catch(function () {
          /* Surface meshing fails on a very small or disordered ligand.
             Inflated spheres are the same idea with no mesh to build. */
          return builder.addRepresentation(component.cell, {
            type: 'spacefill',
            typeParams: {
              alpha: NEON.cloud.alpha, emissive: NEON.cloud.emissive,
              sizeFactor: 1.7, ignoreLight: true
            },
            color: 'uniform',
            colorParams: { value: NEON.neon }
          }).catch(function () { /* no halo; the sticks still glow */ });
        });
      });
    } catch (err) { /* no halo; the sticks still glow */ }
  }

  /* Frame linking: the viewer and the time-series plots share one x axis
   * (simulated time), so moving in either should move the other. Kept
   * deliberately simple -- a vertical line on each time plot, updated as the
   * viewer's frame changes. */

  function timePlots() {
    return drawn.filter(function (node) {
      const layout = node.layout;
      return layout && layout.xaxis && layout.xaxis.title &&
        String(layout.xaxis.title.text || '').indexOf('Time') === 0;
    });
  }

  /* The cursor is appended to whatever shapes a panel already has, not
     substituted for them. Several binding panels carry their own: the shaded
     equilibration head, the 5 A drift gate, the quarter-retained contact
     threshold. Replacing the list wiped all of those the first time the slider
     moved, so a panel lost its own annotations as soon as anybody scrubbed it. */
  const baseShapes = new WeakMap();

  /* Which time panels are actually on screen.
     
     A relayout costs the same whether anyone can see the panel or not, and on a
     reference-state run page there are ten time panels down a page several
     screens long: while someone watches the viewer at the top, nine of them are
     being redrawn out of view. The observer is created once and panels are only
     added to it, so a page with no IntersectionObserver (or none yet observed)
     falls back to redrawing everything, which is what it did before. */
  let onScreen = null;
  function visiblePlots() {
    const plots = timePlots();
    if (typeof IntersectionObserver !== 'function') return plots;
    if (!onScreen) {
      onScreen = new WeakMap();
      const observer = new IntersectionObserver(function (entries) {
        entries.forEach(function (entry) {
          onScreen.set(entry.target, entry.isIntersecting);
        });
      }, { rootMargin: '120px' });
      plots.forEach(function (node) { observer.observe(node); });
      return plots;                /* nothing reported yet: draw them all once */
    }
    const seen = plots.filter(function (node) { return onScreen.has(node); });
    if (!seen.length) return plots;
    return seen.filter(function (node) { return onScreen.get(node); });
  }

  /* The cursor, coalesced to its own rate.
     
     The frame clock steps the model up to twenty-five times a second, and the
     cursor does not need to: a 1.5px line moving across a 600px panel is at the
     same pixel for several frames at a time, and a relayout per panel per frame
     was the second half of the stutter the first half of this file explains.
     Eight a second is beyond what the eye resolves on a line that thin, and it
     is the frequency Plotly can sustain across ten panels without taking the
     main thread away from the compositor.
     
     Coalesced rather than dropped: the last requested position is always drawn,
     so letting go of the slider never leaves the cursor behind the model. */
  const CURSOR_MS = 125;
  let cursorPending = null;
  let cursorTimer = null;

  function drawCursor(ns) {
    visiblePlots().forEach(function (node) {
      if (!baseShapes.has(node)) {
        const existing = (node.layout && node.layout.shapes) || [];
        baseShapes.set(node, existing.map(function (s) { return s; }));
      }
      const shapes = baseShapes.get(node).concat([{
        type: 'line', x0: ns, x1: ns, y0: 0, y1: 1, yref: 'paper',
        line: { color: '#1a1a2e', width: 1.5 }
      }]);
      Plotly.relayout(node, { shapes: shapes });
    });
  }

  function markFrame(ns) {
    cursorPending = ns;
    if (cursorTimer) return;
    cursorTimer = setTimeout(function () {
      cursorTimer = null;
      const value = cursorPending;
      cursorPending = null;
      if (value !== null) drawCursor(value);
    }, CURSOR_MS);
  }

  /* ------------------------------------------------------- frame slider */

  /* Step a multi-model structure WITHOUT reloading it.

     The frame is a parameter on whichever state cell carries a modelIndex, found
     by that rather than by transformer name: the name differs between the
     trajectory and the model-from-trajectory transforms, and this only needs the
     thing with the knob on it. Updating it costs a parameter commit and never
     touches the camera, so zoom and orientation stay where the user left them.
     Loading a frame as a structure instead would make Mol* tear the scene down
     and re-frame it, which is what makes a scrub flash and jump. */
  function setFrame(viewer, index) {
    try {
      const plugin = viewer.plugin;
      if (!plugin || !plugin.state || !plugin.build) return false;
      const cells = [];
      plugin.state.data.cells.forEach(function (cell) {
        const params = cell.transform && cell.transform.params;
        if (params && typeof params.modelIndex === 'number') cells.push(cell);
      });
      if (!cells.length) return false;
      const update = plugin.build();
      cells.forEach(function (cell) {
        update.to(cell.transform.ref).update(function (old) {
          old.modelIndex = index;
        });
      });
      update.commit();
      return true;
    } catch (err) {
      return false;
    }
  }

  /* Wire the slider, the readout, the play button and the plot cursor together.

     Returns true when the slider is driving playback, so the caller knows not
     to also start Mol*'s own animation loop: two things stepping the same model
     index fight, and the slider's position then disagrees with what is drawn. */
  function wireSlider(viewer, frames, duration) {
    const slider = document.getElementById('scrub');
    if (!slider || frames < 2) return false;
    const readout = document.getElementById('scrub-readout');
    const play = document.getElementById('scrub-play');

    function show() {
      const index = parseInt(slider.value, 10);
      const ns = duration > 0 ? (duration * index) / (frames - 1) : index;
      if (readout) {
        readout.textContent = duration > 0
          ? ns.toFixed(duration < 10 ? 3 : 1) + ' ns'
          : 'frame ' + index;
      }
      markFrame(ns);
      /* No debounce: a model switch is a parameter update rather than a load,
         so dragging can drive it directly. */
      setFrame(viewer, index);
    }
    slider.addEventListener('input', show);

    /* Set only by the person, never by the page. Without it, anything that
       restarts playback would override a deliberate pause. */
    let paused = window.matchMedia
      && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

    /* Driven by requestAnimationFrame with a time accumulator rather than by
       setInterval, for three reasons that all showed up as stutter. A timer
       fires on its own schedule and the browser composites on its own, so every
       few frames a model update lands just after a paint and is seen a frame
       late. A timer also keeps firing in a backgrounded tab, where it is
       throttled to once a second and then delivers a burst of queued work on
       return. And a fractional interval cannot be expressed at all: 500 frames
       over 20 s is 40 ms, which is fine, but 180 frames over 8 s is 44.4 ms and
       setInterval rounds it, so the pass drifts.
       
       The accumulator carries the remainder between callbacks, so the playback
       runs at the real rate on any display and skips rather than stalls if a
       callback is late. */
    let frameClock = null;
    let lastTick = 0;
    let carry = 0;

    function stop(byUser) {
      if (byUser) paused = true;
      if (frameClock) { cancelAnimationFrame(frameClock); frameClock = null; }
      lastTick = 0;
      carry = 0;
      if (play) {
        play.textContent = 'Play';
        play.setAttribute('aria-label', 'Play the trajectory');
      }
    }
    function start(byUser) {
      if (byUser) paused = false;
      if (frameClock || paused) return;
      if (play) {
        play.textContent = 'Pause';
        play.setAttribute('aria-label', 'Pause the trajectory');
      }
      const perMs = frames / (playbackSeconds(frames) * 1000);
      lastTick = 0;
      carry = 0;
      frameClock = requestAnimationFrame(function tick(now) {
        if (!frameClock) return;
        const delta = lastTick ? Math.min(now - lastTick, 250) : 0;
        lastTick = now;
        carry += delta * perMs;
        const advance = Math.floor(carry);
        if (advance >= 1) {
          carry -= advance;
          slider.value = String((parseInt(slider.value, 10) + advance) % frames);
          show();
        }
        frameClock = requestAnimationFrame(tick);
      });
    }
    if (play) {
      play.addEventListener('click', function () {
        if (frameClock) { stop(true); } else { start(true); }
      });
    }
    /* Dragging pauses: the person has taken the wheel, and a loop that keeps
       advancing underneath them is fighting the hand on the control. */
    slider.addEventListener('pointerdown', function () { stop(true); });

    /* Published for the structure toggles. `pause` reports whether it was
       playing, so a caller can leave it paused if the person already had. */
    playback = {
      pause: function () {
        const wasPlaying = Boolean(frameClock);
        if (wasPlaying) stop(false);
        return wasPlaying;
      },
      resume: function (wasPlaying) { if (wasPlaying) start(false); },
    };

    show();
    start(false);
    return true;
  }

  /* Start the trajectory looping as soon as it is loaded.

     Mol* registers its animations by name and does not export the model-index
     one on the global bundle, so it is looked up rather than imported -- which
     also means a version that renames or drops it degrades to a static
     structure with working manual controls, instead of throwing inside the
     load chain and losing the viewer entirely.

     Skipped for anyone who has asked their system for reduced motion. A
     3D scene that starts moving on its own is exactly what that setting is
     about, and the playback controls are still there for them. */
  function autoplay(viewer, frames) {
    if (window.matchMedia
        && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
      return false;
    }
    try {
      const manager = viewer.plugin.managers.animation;
      const animation = manager.animations.find(function (a) {
        return a.name === 'built-in.animate-model-index';
      });
      if (!animation) return false;
      manager.play(animation, {
        mode: { name: 'loop', params: { direction: 'forward' } },
        duration: { name: 'fixed',
                    params: { durationInS: playbackSeconds(frames) } }
      });
      return true;
    } catch (err) {
      return false;      /* a static structure is a fine outcome; a broken page is not */
    }
  }

  /* Mol* injects its own logo as a link with no text and no label, which
   * axe-core reports as a serious link-name violation. It is the library's
   * markup rather than ours, so it is labelled after the viewer builds rather
   * than fixed at the source. Named rather than hidden: it is a real link to
   * the project and crediting it is the point of it being there. */
  function nameTheMolstarLogo(host) {
    /* Named in the same task the logo is inserted in, rather than on the next
     * poll. Mol* builds its viewport controls during creation, and with the
     * viewer now started after first paint, axe-core sampled the page inside
     * the window between the logo appearing and a 100 ms poll finding it: the
     * check reported the exact violation this function exists to prevent. The
     * observer closes that window and stops the moment it has named the link,
     * so it is still not left running for the life of the page; the bounded
     * poll stays as the fallback. */
    let done = false;
    let observer = null;

    function label() {
      const logo = host.querySelector('.msp-logo');
      if (!logo) return false;
      if (!logo.getAttribute('aria-label')) {
        logo.setAttribute('aria-label', 'Mol* molecular viewer, opens molstar.org');
      }
      if (!logo.getAttribute('title')) {
        logo.setAttribute('title', 'Mol* molecular viewer');
      }
      return true;
    }

    function stop() {
      done = true;
      if (observer) { observer.disconnect(); observer = null; }
    }

    if (label()) return;
    if (typeof MutationObserver === 'function') {
      observer = new MutationObserver(function () { if (label()) stop(); });
      observer.observe(host, {childList: true, subtree: true});
    }
    let attempts = 0;
    (function attempt() {
      if (done) return;
      if (label() || ++attempts >= 40) { stop(); return; }
      setTimeout(attempt, 100);
    })();
  }

  /* Section 6.2's membrane toggle. The bilayer is a second structure loaded
   * beside the receptor rather than part of it: the viewer bundle holds
   * protein and ligand only, and merging 476 phosphorus atoms into the
   * topology every measurement is taken against is not worth the convenience.
   *
   * Loaded on demand and deleted again, rather than loaded and hidden. The
   * files are small, so there is nothing to save by keeping them resident, and
   * deleting the subtree by the ref `loadTrajectory` returns needs no
   * visibility API, which is the part of Mol*'s surface most likely to move
   * between versions. */
  /* Section 6.2's reference overlay, one ghost per state.
   *
   * Ghosted means a translucent cartoon in that state's colour, so the eye
   * reads it as background against the run's own motif colours. The colours
   * are the ones the badges and the state map use, which is what lets a reader
   * carry "green is active" from one part of the page to the other.
   *
   * The structure is static: it was fitted on the core of the frame this
   * viewer opens on, by `tools/analyse_run.py`, and checked there against the
   * motif distances the analysis stored. Drift away from it as the slider runs
   * is the run moving, which is the question the page is about.
   */
  function overlayColour(state) {
    const root = getComputedStyle(document.documentElement);
    const value = cssColourToInt(root.getPropertyValue('--state-' + state));
    return value === null ? 0x8a94a6 : value;
  }

  /* Every component ref currently in the hierarchy. Taken before a load and
     differenced after it, which is how the newly added structure's components
     are found without depending on the shape of what `loadStructureFromUrl`
     returns. */
  function componentRefs(plugin) {
    const refs = [];
    try {
      (plugin.managers.structure.hierarchy.current.structures || [])
        .forEach(function (ref) {
          (ref.components || []).forEach(function (component) {
            if (component.cell && component.cell.transform) {
              refs.push(component.cell.transform.ref);
            }
          });
        });
    } catch (err) { /* nothing loaded yet */ }
    return refs;
  }

  function ghostComponents(plugin, before, colour) {
    const seen = {};
    before.forEach(function (ref) { seen[ref] = true; });
    let touched = 0;
    try {
      const update = plugin.build();
      (plugin.managers.structure.hierarchy.current.structures || [])
        .forEach(function (ref) {
          (ref.components || []).forEach(function (component) {
            const key = component.cell && component.cell.transform
              && component.cell.transform.ref;
            if (!key || seen[key]) return;
            (component.representations || []).forEach(function (rep) {
              if (!rep.cell) return;
              update.to(rep.cell).update(function (old) {
                if (!old.type) return;
                old.type.params = old.type.params || {};
                /* Translucent, so the run reads in front of it. */
                old.type.params.alpha = 0.3;
                old.colorTheme = { name: 'uniform',
                                   params: { value: colour } };
                touched += 1;
              });
            });
          });
        });
      update.commit();
    } catch (err) {
      return 0;
    }
    return touched;
  }

  function wireOverlayToggles(viewer, say) {
    const toggles = document.querySelectorAll('.overlay-toggle');
    if (!toggles.length) return;

    toggles.forEach(function (toggle) {
      let loaded = null;
      let busy = false;
      const state = toggle.dataset.state || 'active';

      toggle.addEventListener('change', async function () {
        if (busy) return;
        busy = true;
        toggle.disabled = true;
        /* Paused for the load and resumed after it, so the structure arrives
           in a second rather than in ten. */
        const wasPlaying = playback.pause();
        if (say) {
          say(toggle.checked ? 'Loading the ' + state + ' reference\u2026'
                             : 'Removing the overlay\u2026');
        }
        try {
          if (toggle.checked) {
            /* Coloured after the load, not through `representationParams`.
             * This bundle ignores that option: the ghosts came back in Mol*'s
             * default chain colours at full opacity, which is two bright
             * proteins on top of the run rather than a ghost behind it.
             * Differencing the component refs across the load is what finds
             * the new structure's representations to retheme. */
            const before = componentRefs(viewer.plugin);
            loaded = await viewer.loadStructureFromUrl(
              toggle.dataset.structure, 'pdb', false);
            ghostComponents(viewer.plugin, before, overlayColour(state));
            toggle.dataset.loaded = 'yes';
            if (say) say('Ghosted the nearest ' + state
                         + ' reference, fitted on the core of the first frame.');
          } else if (loaded && loaded.model) {
            await viewer.plugin.state.data.build()
              .delete(loaded.model.ref).commit();
            loaded = null;
            toggle.dataset.loaded = 'no';
            if (say) say('Overlay hidden.');
          }
        } catch (err) {
          toggle.checked = false;
          toggle.disabled = true;
          if (say) say('The overlay could not be drawn: ' + err, true);
          return;
        } finally {
          busy = false;
          playback.resume(wasPlaying);
        }
        toggle.disabled = false;
      });
    });
  }

  /* Point the camera along the bilayer rather than down through it.
   *
   * A membrane run opens on whatever view Mol* chose when it framed the
   * structure, which is as likely to be straight down the bilayer normal as
   * anything else: the membrane is then a disc the receptor is buried in and
   * the thing a reader is being shown cannot be seen at all. Edge on, it is a
   * slab with the receptor through it.
   *
   * Both vectors are measured by `tools/pack_run_fxa.py` and travel in the
   * results file, because neither is something a page can guess. The normal is
   * the axis the bilayer's phosphorus atoms are least spread along, and which
   * end of it is extracellular comes from helix 8, which lies on the
   * intracellular face: `up` points at the extracellular side, so the
   * extracellular loops are at the top of the viewport.
   *
   * Everything here is optional and guarded. A page with no orientation, or a
   * Mol* whose camera does not expose a snapshot, keeps the default view: a
   * viewer pointing somewhere arbitrary is a worse outcome than one pointing
   * where Mol* put it. */
  /* A copy, not the camera's own object.
     
     `getSnapshot` hands back state that Mol* goes on mutating, so a reference
     kept across a load is not a record of the view before it: the first
     version of this restored the framing to exactly what it was trying to
     undo, and the receptor stayed a speck in the middle of a bilayer four
     times its width. The three vectors are copied element by element because
     they are typed arrays. */
  function cameraSnapshot(viewer) {
    try {
      const camera = viewer.plugin && viewer.plugin.canvas3d
        && viewer.plugin.canvas3d.camera;
      if (!camera || typeof camera.getSnapshot !== 'function') return null;
      const live = camera.getSnapshot();
      const copy = Object.assign({}, live);
      ['position', 'target', 'up'].forEach(function (key) {
        if (live[key] && live[key].length >= 3) {
          copy[key] = [live[key][0], live[key][1], live[key][2]];
        }
      });
      return copy;
    } catch (err) {
      return null;
    }
  }

  function restoreCamera(viewer, snapshot) {
    if (!snapshot) return false;
    try {
      viewer.plugin.canvas3d.camera.setState(snapshot);
      return true;
    } catch (err) {
      return false;
    }
  }

  /* Applied now and again shortly after.
     
     Mol* frames the camera on what it holds when a load resolves, and that
     framing runs on its own schedule *after* the promise settles: a camera set
     inside the continuation is overwritten a moment later, which is why the
     first version of this produced a perfect top-down view of a bilayer that
     was supposed to be edge on. The second application lands after Mol* has
     finished and is what the reader sees. Both are guarded, so a bundle whose
     camera does not expose a snapshot simply keeps its own view. */
  function orientToMembraneSoon(viewer, host) {
    /* Applied, then checked and applied again if it did not hold.
       
       Mol* frames the camera on what it holds when a load resolves, and it
       animates into that framing over the next few hundred milliseconds: a
       camera set inside the continuation is interpolated away a moment later,
       which is why the first two attempts at this produced a perfect top-down
       view of a bilayer that was supposed to be edge on. Checking the result
       and reapplying is what makes it stick without guessing a delay long
       enough to outlast an animation whose length is Mol*'s business. Three
       looks over two seconds, and it stops as soon as the view is the one
       asked for, so a reader who takes hold of the camera is not fighting it. */
    let plan;
    try {
      plan = JSON.parse(host && host.dataset.orientation || 'null');
    } catch (err) {
      return;
    }
    if (!plan || !plan.up) return;
    orientToMembrane(viewer, host);
    let attempts = 0;
    const settle = setInterval(function () {
      attempts += 1;
      const snapshot = cameraSnapshot(viewer);
      const held = snapshot && snapshot.up
        && Math.abs(snapshot.up[0] - plan.up[0]) < 0.01
        && Math.abs(snapshot.up[1] - plan.up[1]) < 0.01
        && Math.abs(snapshot.up[2] - plan.up[2]) < 0.01;
      if (held || attempts >= 3) { clearInterval(settle); return; }
      orientToMembrane(viewer, host);
    }, 700);
  }

  function orientToMembrane(viewer, host) {
    if (!host || !host.dataset.orientation) return false;
    let plan;
    try {
      plan = JSON.parse(host.dataset.orientation);
    } catch (err) {
      return false;
    }
    if (!plan || !plan.up || !plan.along) return false;
    try {
      const camera = viewer.plugin && viewer.plugin.canvas3d
        && viewer.plugin.canvas3d.camera;
      if (!camera || typeof camera.getSnapshot !== 'function') return false;
      const snapshot = camera.getSnapshot();
      const target = snapshot.target;
      /* The distance Mol* already chose, so the structure stays framed: this
         moves the camera around the target rather than towards it. */
      const dx = snapshot.position[0] - target[0];
      const dy = snapshot.position[1] - target[1];
      const dz = snapshot.position[2] - target[2];
      const radius = Math.sqrt(dx * dx + dy * dy + dz * dz)
        || (snapshot.radius || 50) * 3;
      camera.setState(Object.assign({}, snapshot, {
        position: [target[0] + plan.along[0] * radius,
                   target[1] + plan.along[1] * radius,
                   target[2] + plan.along[2] * radius],
        up: [plan.up[0], plan.up[1], plan.up[2]],
      }));
      return true;
    } catch (err) {
      return false;
    }
  }

  /* The bilayer, drawn so the receptor is still visible through it.
   *
   * Mol*'s default preset for a structure of 476 lipids and no polymer is
   * spheres, and 24,752 of them is a wall: the first version drew a perfect
   * edge-on bilayer with the receptor completely hidden inside it. Lines show
   * the chains, which is what a membrane looks like, and cost almost nothing
   * to draw. Translucent on top of that, so the cartoon reads in front.
   *
   * Differenced across the load like the ghosts, for the same reason: it is how
   * the new structure's representations are found without depending on the
   * shape of what the loader returns. Any failure leaves the default, which is
   * ugly rather than broken.
   */
  function thinComponents(plugin, before) {
    const seen = {};
    before.forEach(function (ref) { seen[ref] = true; });
    let touched = 0;
    try {
      const update = plugin.build();
      (plugin.managers.structure.hierarchy.current.structures || [])
        .forEach(function (ref) {
          (ref.components || []).forEach(function (component) {
            const key = component.cell && component.cell.transform
              && component.cell.transform.ref;
            if (!key || seen[key]) return;
            (component.representations || []).forEach(function (rep) {
              if (!rep.cell) return;
              update.to(rep.cell).update(function (old) {
                if (!old.type) return;
                old.type.name = 'line';
                old.type.params = old.type.params || {};
                /* Transparent, so the receptor reads through the membrane
                   rather than behind it. The bilayer is context: a reader is
                   looking at the protein in it, not at the lipids. */
                old.type.params.alpha = 0.3;
                touched += 1;
              });
            });
          });
        });
      update.commit();
    } catch (err) {
      return 0;
    }
    return touched;
  }

  function wireLipidToggle(viewer, say, host) {
    const toggle = document.getElementById('lipid-toggle');
    if (!toggle) return;
    let loaded = null;
    let busy = false;

    toggle.addEventListener('change', async function () {
      if (busy) return;
      busy = true;
      toggle.disabled = true;
      /* Paused for the load, for the reason the overlay toggle gives: the
         bilayer is a second trajectory and loading it while the first one
         plays takes the thread from both. */
      const wasPlaying = playback.pause();
      if (say) {
        say(toggle.checked ? 'Loading the bilayer\u2026' : 'Hiding the bilayer\u2026');
      }
      try {
        if (toggle.checked) {
          /* A trajectory when there is one, a structure when there is not.
             The bilayer is the whole lipid now rather than its phosphorus,
             which is 24,752 atoms: two megabytes for one frame and forty-eight
             for five hundred, so it is drawn at the first frame and the
             control says so. Without a trajectory there is nothing with a
             model index in it, which also keeps it out of the scrubber's way. */
          /* Taken before the load and restored after it. Mol* frames the
             camera on everything it holds, so adding a bilayer four times the
             width of the receptor zooms out until the receptor is a speck in
             the middle of it. The run is what the page is about. */
          const framing = cameraSnapshot(viewer);
          const before = componentRefs(viewer.plugin);
          loaded = toggle.dataset.trajectory
            ? await viewer.loadTrajectory({
                model: { kind: 'model-url', url: toggle.dataset.structure,
                         format: 'pdb' },
                coordinates: { kind: 'coordinates-url',
                               url: toggle.dataset.trajectory,
                               format: 'xtc', isBinary: true },
                preset: 'default'
              })
            : await viewer.loadStructureFromUrl(
                toggle.dataset.structure, 'pdb', false);
          /* Said out loud. A control that does something invisible from this
           * camera angle, which the bilayer is until you rotate, reads as a
           * control that did nothing. The attribute is the same fact in a
           * form a test can assert on without reaching into Mol*. */
          thinComponents(viewer.plugin, before);
          restoreCamera(viewer, framing);
          toggle.dataset.loaded = 'yes';
          orientToMembraneSoon(viewer, host);
          if (say) {
            say(toggle.dataset.trajectory
              ? 'Bilayer drawn around the receptor, moving with it.'
              : 'Bilayer drawn around the receptor, at the first frame.');
          }
        } else if (loaded && loaded.model) {
          await viewer.plugin.state.data.build()
            .delete(loaded.model.ref).commit();
          loaded = null;
          toggle.dataset.loaded = 'no';
          if (say) say('Bilayer hidden.');
        }
      } catch (err) {
        /* A toggle that silently does nothing is worse than one that says it
         * could not. The control is left off and disabled rather than
         * inviting a second click that will fail the same way. */
        toggle.checked = false;
        toggle.disabled = true;
        if (say) say('The bilayer could not be drawn: ' + err, true);
        return;
      } finally {
        busy = false;
        playback.resume(wasPlaying);
        if (!toggle.disabled) toggle.disabled = false;
      }
      toggle.disabled = false;
    });
  }

  function initViewer() {
    const host = document.getElementById('molstar-viewer');
    const status = document.getElementById('viewer-status');
    if (!host || typeof molstar === 'undefined') return;

    /* Started before the viewer is created, not after it resolves: the logo
     * goes in during creation, so a labeller that waits for the promise has
     * already lost a window where the link has no name. */
    nameTheMolstarLogo(host);

    const structureUrl = host.dataset.structure;
    const trajectoryUrl = host.dataset.trajectory;
    const duration = parseFloat(host.dataset.duration || '0');
    const frames = parseInt(host.dataset.frames || '0', 10);

    function say(message, isError) {
      if (!status) return;
      status.textContent = message;
      status.classList.toggle('is-error', Boolean(isError));
    }

    molstar.Viewer.create(host, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: true,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      viewportShowAnimation: true,
      pdbProvider: 'rcsb',
      emdbProvider: 'rcsb'
    }).then(function (viewer) {
      /* Labelled as soon as the viewer exists, not after the structure loads.
       * Mol* builds its viewport controls, logo included, during creation, so
       * a structure that fails to load still leaves the logo on the page and
       * the labelling inside the load handler never ran: axe reported the
       * link-name violation again on exactly that path. */
      nameTheMolstarLogo(host);
      /* Topology plus coordinates: Mol* reads the PDB for connectivity and the
       * XTC for the frames, which is exactly how the .fxa stores them. */
      return viewer.loadTrajectory({
        model: { kind: 'model-url', url: structureUrl, format: 'pdb' },
        coordinates: { kind: 'coordinates-url', url: trajectoryUrl,
                       format: 'xtc', isBinary: true },
        preset: 'default'
      }).then(function () {
        /* The slider is the one thing stepping the model, where it exists.
         *
         * It replaces both halves of what was here: Mol*'s own model-index
         * animation, and a 250 ms poll that watched the frame change so the plot
         * cursor could follow. Two things driving one model index fight, and the
         * slider's position then disagrees with what is drawn; and a poll could
         * only ever read the frame, never set it, so there was no way to go to a
         * frame a plot made you curious about. Driving it from one control makes
         * the structure and the cursor the same fact.
         *
         * Mol*'s animation stays as the fallback for a page without the slider,
         * which is any single-frame structure. */
        applyNeon(viewer);
        colourByMotif(viewer, host);
        wireLipidToggle(viewer, say, host);
        wireOverlayToggles(viewer, say);
        /* The bilayer is on by default on a membrane run, so the control that
           is already checked has to be acted on rather than waited for: a
           `change` listener never fires for a box that arrived checked.
           
           Deferred to a later task rather than dispatched here. Fired inline
           it starts a second structure load from inside the continuation of
           the first one, and Mol* never resolves it: the toggle sat disabled
           and the bilayer never appeared, with no error anywhere, while the
           same control clicked by hand a second later worked every time. */
        const lipids = document.getElementById('lipid-toggle');
        if (lipids && lipids.checked) {
          setTimeout(function () {
            lipids.dispatchEvent(new Event('change', { bubbles: true }));
          }, 0);
        } else {
          orientToMembraneSoon(viewer, host);
        }
        const scrubbing = wireSlider(viewer, frames, duration);
        const playing = scrubbing || (frames > 1 && autoplay(viewer, frames));
        const pass = Math.round(playbackSeconds(frames));
        say(frames
          ? frames + ' frames'
              + (scrubbing
                  ? ', on a ' + pass + ' second loop. Drag the slider to a frame and every time plot marks it.'
                  : (playing
                      ? ', playing on a ' + pass + ' second loop. The controls under the viewport pause it.'
                      : ' loaded. Use the playback controls under the viewport.'))
          : 'Structure loaded.');
      });
    }).catch(function (err) {
      say('The 3D viewer could not load this structure (' + err + '). '
        + 'The plots below are unaffected.', true);
    });
  }

  /* ----------------------------------------------------------- boot */

  /* The viewer waits for the first paint rather than racing it.
   *
   * `molstar.Viewer.create` builds a WebGL context and mounts its whole UI on
   * the main thread, so calling it straight out of DOMContentLoaded holds the
   * compositor until that finishes: measured in headless Chromium, every asset
   * had arrived by 50 ms and DOMContentLoaded fired at about 0.5 s, while first
   * contentful paint landed at 3.2 s on the runs where the context fell back to
   * software. The page is complete and readable before any of that: section
   * 6.3's budget is about when a reader sees something, and the plots and the
   * viewer are explicitly allowed to arrive after it.
   *
   * Two frames rather than one. The first callback runs before the paint it was
   * scheduled against is composited, so work started there is still in front of
   * it; the second runs after. The timeout is the fallback for a backgrounded
   * tab, where rAF never fires and the viewer would otherwise never start. */
  function afterFirstPaint(run) {
    if (typeof requestAnimationFrame !== 'function') { setTimeout(run, 0); return; }
    let started = false;
    const once = function () { if (!started) { started = true; run(); } };
    requestAnimationFrame(function () { requestAnimationFrame(once); });
    setTimeout(once, 1000);
  }

  document.addEventListener('DOMContentLoaded', function () {
    initDropZone();
    /* Both of the heavy ones wait for the paint now. The viewer has since
       4bp; the plots joined it when a page grew to thirty-five of them. */
    afterFirstPaint(initPlots);
    afterFirstPaint(initViewer);
  });
})();

/* ------------------------------------------------------- re-analysis -----
 * Fire a bounded server-side job, then poll. The work happens in a detached
 * subprocess, so this is a status poll rather than a long-held request.
 */

(function () {
  'use strict';

  const POLL_MS = 1500;
  const PLOT_CONFIG = {
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d']
  };

  function init() {
    const card = document.getElementById('reanalyse-card');
    if (!card) return;

    const button = document.getElementById('re-run');
    const status = document.getElementById('re-status');
    const output = document.getElementById('re-results');
    const selection = document.getElementById('re-selection');
    const selectionB = document.getElementById('re-selection-b');
    const selectionBField = document.getElementById('re-selection-b-field');

    function chosen() {
      return Array.from(document.querySelectorAll('input[name="re-metric"]:checked'))
        .map(function (el) { return el.value; });
    }

    /* The second selection only means anything for a distance. */
    function syncDistanceField() {
      selectionBField.hidden = chosen().indexOf('distance') === -1;
    }
    document.querySelectorAll('input[name="re-metric"]').forEach(function (el) {
      el.addEventListener('change', syncDistanceField);
    });
    syncDistanceField();

    function say(message, kind) {
      status.textContent = message;
      status.className = 'md-reanalyse-status' + (kind ? ' is-' + kind : '');
    }

    function draw(panels) {
      output.innerHTML = '';
      if (!panels || !panels.length) {
        say('That produced nothing to plot.', 'error');
        return;
      }
      panels.forEach(function (panel) {
        const wrap = document.createElement('div');
        wrap.className = 'md-panel md-reanalyse-panel';
        const heading = document.createElement('h4');
        heading.textContent = panel.title;
        const target = document.createElement('div');
        target.className = 'md-plot';
        target.id = 'plot-' + panel.id;
        wrap.appendChild(heading);
        wrap.appendChild(target);
        output.appendChild(wrap);
        Plotly.newPlot(target, panel.figure.data, panel.figure.layout, PLOT_CONFIG);
      });
    }

    let timer = null;

    function poll() {
      fetch(card.dataset.status)
        .then(function (r) { return r.json(); })
        .then(function (data) {
          if (data.status === 'running') {
            say('Working' + (data.seconds ? ' (' + data.seconds + 's)' : '') + '…');
            timer = setTimeout(poll, POLL_MS);
            return;
          }
          button.disabled = false;
          if (data.status === 'ready') {
            const m = data.metrics || {};
            say('Done — ' + m.n_atoms + ' atoms, ' + m.n_frames + ' frames.', 'ok');
            draw(data.figures);
          } else {
            say(data.message || 'That did not work.', 'error');
          }
        })
        .catch(function () {
          button.disabled = false;
          say('Lost contact with the server.', 'error');
        });
    }

    button.addEventListener('click', function () {
      const metrics = chosen();
      if (!metrics.length) {
        say('Choose at least one metric.', 'error');
        return;
      }
      clearTimeout(timer);
      button.disabled = true;
      say('Starting…');
      output.innerHTML = '';

      fetch(card.dataset.endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          metrics: metrics,
          selection: selection.value,
          selection_b: selectionB ? selectionB.value : ''
        })
      })
        .then(function (r) { return r.json().then(function (d) { return [r.status, d]; }); })
        .then(function (pair) {
          const code = pair[0], data = pair[1];
          if (code === 202) { timer = setTimeout(poll, POLL_MS); return; }
          button.disabled = false;
          say(data.message || 'That request was refused.', 'error');
        })
        .catch(function () {
          button.disabled = false;
          say('Could not reach the server.', 'error');
        });
    });
  }

  document.addEventListener('DOMContentLoaded', init);
})();

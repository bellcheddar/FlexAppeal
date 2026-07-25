/* FlexAppeal front end.
 *
 * Vanilla JS, no framework, no build step. It does exactly two things:
 *
 *   1. switches the visible input on the structure-source picker
 *   2. keeps the live readout panel and field visibility in sync with the form
 *
 * Note what it deliberately does NOT do: decide which fields apply. That is
 * decided by the `requires` predicates in options.py, evaluated server-side, and
 * returned as an `active` list by /api/estimate. Reimplementing that logic here
 * would mean two predicate engines that agree right up until they do not.
 */

(function () {
  'use strict';

  const DEBOUNCE_MS = 300;

  /* ------------------------------------------------------------ helpers */

  function debounce(fn, ms) {
    let timer = null;
    return function () {
      clearTimeout(timer);
      timer = setTimeout(fn, ms);
    };
  }

  function formatInt(n) {
    return Number(n).toLocaleString('en-GB');
  }

  /* --------------------------------------------- structure source picker */

  function initSourcePicker() {
    const radios = document.querySelectorAll('input[name="input_source"][type="radio"]');
    const panels = document.querySelectorAll('.md-source-inputs [data-source]');
    if (!radios.length || !panels.length) return;

    function sync() {
      const checked = document.querySelector('input[name="input_source"]:checked');
      const value = checked ? checked.value : 'upload';
      panels.forEach(function (panel) {
        const sources = panel.dataset.source.split(/\s+/);
        panel.hidden = sources.indexOf(value) === -1;
      });
    }

    radios.forEach(function (radio) { radio.addEventListener('change', sync); });
    sync();
  }

  /* ------------------------------------------------- form serialisation */

  /* Distinguishing a multiselect from a boolean by DOM position rather than by
   * a duplicated widget table: a multiselect renders its options inside
   * .md-checkgroup, a boolean renders a single .md-switch. */
  function serialise(form) {
    const data = {};

    form.querySelectorAll('[name]').forEach(function (el) {
      const name = el.name;
      if (!name || name.startsWith('_')) return;

      if (el.type === 'checkbox') {
        if (el.closest('.md-checkgroup')) {
          if (!Array.isArray(data[name])) data[name] = [];
          if (el.checked) data[name].push(el.value);
        } else {
          data[name] = el.checked;
        }
      } else if (el.type === 'radio') {
        if (el.checked) data[name] = el.value;
      } else if (el.type === 'file') {
        /* skip -- the file is already on the server for this session */
      } else {
        data[name] = el.value;
      }
    });

    return data;
  }

  /* ------------------------------------------------------- live readouts */

  function initReadout() {
    const form = document.getElementById('prepare-form');
    const readout = document.getElementById('readout');
    if (!form || !readout) return;

    const cells = {};
    readout.querySelectorAll('[data-readout]').forEach(function (el) {
      cells[el.dataset.readout] = el;
    });

    let inFlight = null;

    function applyActive(active) {
      const set = new Set(active);
      document.querySelectorAll('.md-field[data-option]').forEach(function (field) {
        /* A field with no `requires` is always active and is never in the list's
         * way; only hide ones the server actually excluded. */
        if (!field.dataset.requires) return;
        field.hidden = !set.has(field.dataset.option);
      });
    }

    function applyIssues(errors, warnings) {
      document.querySelectorAll('.md-field-error').forEach(function (el) {
        el.hidden = true;
        el.textContent = '';
        const field = el.closest('.md-field');
        if (field) field.classList.remove('has-error');
      });

      errors.concat(warnings).forEach(function (issue) {
        if (!issue.option) return;
        const el = document.querySelector('[data-error-for="' + issue.option + '"]');
        if (!el) return;
        el.textContent = issue.message;
        el.hidden = false;
        const field = el.closest('.md-field');
        if (field && errors.indexOf(issue) !== -1) field.classList.add('has-error');
      });

      const panel = cells.issues;
      if (!panel) return;
      panel.innerHTML = '';
      errors.forEach(function (issue) {
        const div = document.createElement('div');
        div.className = 'md-readout-issue is-error';
        div.textContent = issue.message;
        panel.appendChild(div);
      });
      warnings.slice(0, 3).forEach(function (issue) {
        const div = document.createElement('div');
        div.className = 'md-readout-issue is-warning';
        div.textContent = issue.message;
        panel.appendChild(div);
      });
    }

    function refresh() {
      const payload = serialise(form);
      payload._estimated_atoms = Number(readout.dataset.estimatedAtoms || 0);
      payload._solute_atoms = Number(readout.dataset.soluteAtoms || 0);

      if (inFlight) inFlight.abort();
      const controller = new AbortController();
      inFlight = controller;

      fetch('/api/estimate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
        signal: controller.signal
      })
        .then(function (r) {
          if (!r.ok) throw new Error('estimate failed: ' + r.status);
          return r.json();
        })
        .then(function (data) {
          const d = data.derived;
          const w = data.wall;

          if (cells.ns) cells.ns.textContent = Math.round(d.total_ns) + ' ns';
          if (cells.steps) cells.steps.textContent = formatInt(d.total_steps);
          if (cells.frames) cells.frames.textContent = formatInt(d.traj_frames);
          if (cells.size) cells.size.textContent = d.traj_size_human;
          if (cells.wall) cells.wall.textContent = w.human;
          if (cells.basis) {
            cells.basis.textContent = w.basis === 'estimated'
              ? 'Rough estimate — the bundle benchmarks your hardware before it starts and replaces this with a real figure.'
              : 'Load a structure for a size-aware estimate.';
          }

          applyActive(data.active);
          applyIssues(data.errors, data.warnings);
          readout.classList.remove('is-stale');
        })
        .catch(function (err) {
          if (err.name === 'AbortError') return;
          /* A failed estimate must never block the form -- the server
           * revalidates on submit regardless, so stale numbers are marked
           * rather than hidden. */
          readout.classList.add('is-stale');
        });
    }

    const scheduled = debounce(refresh, DEBOUNCE_MS);
    form.addEventListener('input', scheduled);
    form.addEventListener('change', scheduled);
    refresh();
  }

  /* ----------------------------------------------------------- boot */

  document.addEventListener('DOMContentLoaded', function () {
    initSourcePicker();
    initReadout();
  });
})();

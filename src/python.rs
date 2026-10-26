//! The Python boundary owns bytes; all geometry work executes without the GIL.
use crate::{
    Boundary, coverage_simplify,
    topology::{self, Kind, Shape},
    wkb,
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};
use rayon::prelude::*;
use std::fmt;

/// Index paths from a geometry to its standalone `LinearRing` members.
type RingPaths = Vec<Vec<usize>>;
type Simplified<'py> = (Bound<'py, PyBytes>, RingPaths);

fn value_error(e: impl fmt::Display) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Runs `work` on a dedicated pool of `threads`, or on the global pool.
fn run_in_pool<T: Send>(
    threads: Option<usize>,
    work: impl FnOnce() -> PyResult<T> + Send,
) -> PyResult<T> {
    match threads {
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(value_error)?
            .install(work),
        None => work(),
    }
}

/// Re-tags the lines at `paths` as rings, which WKB cannot express.
fn mark_rings(g: &mut Shape, paths: &[Vec<usize>]) -> PyResult<()> {
    for path in paths {
        let mut part = &mut *g;
        for &i in path {
            part = part
                .parts
                .get_mut(i)
                .ok_or_else(|| value_error("invalid ring path"))?;
        }
        if part.kind != Kind::LineString {
            return Err(value_error("ring path must reference a line"));
        }
        part.kind = Kind::LinearRing;
    }
    Ok(())
}

/// Paths to rings reachable through collections; polygon rings are not reported.
fn ring_paths(g: &Shape, path: &mut Vec<usize>, rings: &mut RingPaths) {
    if g.kind == Kind::LinearRing {
        rings.push(path.clone());
    }
    if g.kind.is_collection() {
        for (i, part) in g.parts.iter().enumerate() {
            path.push(i);
            ring_paths(part, path, rings);
            path.pop();
        }
    }
}

fn simplify_one(
    bytes: &[u8],
    tolerance: f64,
    paths: &[Vec<usize>],
) -> PyResult<(Vec<u8>, RingPaths)> {
    if tolerance.is_nan() || tolerance < 0.0 {
        return Err(value_error("Tolerance must be non-negative"));
    }
    let mut g = topology::wkb::read(bytes).map_err(value_error)?;
    mark_rings(&mut g, paths)?;
    g.check_indexable().map_err(value_error)?;
    let result = topology::simplify(g, tolerance);
    let mut rings = Vec::new();
    ring_paths(&result, &mut Vec::new(), &mut rings);
    Ok((topology::wkb::write(&result), rings))
}

#[pyfunction]
#[pyo3(signature = (inputs, tolerances, rings, threads=None))]
fn topology_wkb(
    py: Python<'_>,
    inputs: Vec<Vec<u8>>,
    tolerances: Vec<f64>,
    rings: Vec<RingPaths>,
    threads: Option<usize>,
) -> PyResult<Vec<Simplified<'_>>> {
    if inputs.len() != tolerances.len() || inputs.len() != rings.len() || threads == Some(0) {
        return Err(value_error("invalid batch lengths or threads"));
    }
    let result = py.detach(move || {
        run_in_pool(threads, || {
            inputs
                .par_iter()
                .zip(&tolerances)
                .zip(&rings)
                .map(|((bytes, &tolerance), paths)| simplify_one(bytes, tolerance, paths))
                .collect::<PyResult<Vec<_>>>()
        })
    })?;
    Ok(result
        .into_iter()
        .map(|(b, r)| (PyBytes::new(py, &b), r))
        .collect())
}

#[pyfunction]
#[pyo3(signature = (inputs, tolerance, simplify_boundary, threads=None))]
#[expect(
    clippy::fn_params_excessive_bools,
    reason = "mirrors the Python keyword argument"
)]
fn simplify_wkb(
    py: Python<'_>,
    inputs: Vec<Vec<u8>>,
    tolerance: f64,
    simplify_boundary: bool,
    threads: Option<usize>,
) -> PyResult<Vec<Bound<'_, PyBytes>>> {
    if threads == Some(0) {
        return Err(value_error("threads must be positive"));
    }
    let boundary = if simplify_boundary {
        Boundary::Simplify
    } else {
        Boundary::Preserve
    };
    let result = py.detach(move || {
        run_in_pool(threads, || {
            let mut geometries = inputs
                .par_iter()
                .map(|b| wkb::read(b).map_err(value_error))
                .collect::<PyResult<Vec<_>>>()?;
            coverage_simplify(&mut geometries, tolerance, boundary);
            Ok(geometries.par_iter().map(wkb::write).collect::<Vec<_>>())
        })
    })?;
    Ok(result.iter().map(|b| PyBytes::new(py, b)).collect())
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(simplify_wkb, module)?)?;
    module.add_function(wrap_pyfunction!(topology_wkb, module)?)?;
    Ok(())
}

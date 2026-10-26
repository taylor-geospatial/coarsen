//! Byte parity against GEOS 3.13.1 fixtures and projected real parcels.
use coarsen::{Boundary, Geometry, coverage_simplify, validate, wkb};

fn unhex(s: &str) -> Vec<u8> {
    s.as_bytes()
        .chunks_exact(2)
        .map(|v| {
            let digits = std::str::from_utf8(v).expect("hex is ASCII");
            u8::from_str_radix(digits, 16).expect("fixture WKB is hex")
        })
        .collect()
}
fn geometries(v: &serde_json::Value) -> Vec<Geometry> {
    v.as_array()
        .expect("fixture geometries are an array")
        .iter()
        .map(|s| {
            let hex = s.as_str().expect("fixture geometries are hex strings");
            wkb::read(&unhex(hex)).expect("fixture WKB decodes")
        })
        .collect()
}
#[test]
fn geos_3131_golden_fixtures() {
    let cases: serde_json::Value = serde_json::from_str(include_str!("golden.json")).unwrap();
    for c in cases.as_array().unwrap() {
        let input = geometries(&c["input"]);
        let mask = validate::invalid_mask(&input);
        let expected: Vec<bool> = c["bad"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_bool().unwrap())
            .collect();
        assert_eq!(mask, expected, "validator {}", c["name"]);
        if mask.iter().any(|&b| b) {
            continue;
        }
        for threads in [1, 8] {
            let mut actual = input.clone();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let tolerance = c["tolerance"].as_f64().unwrap();
            pool.install(|| coverage_simplify(&mut actual, tolerance, Boundary::Simplify));
            let expected = geometries(&c["output"]);
            for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
                assert_eq!(
                    wkb::write(a),
                    wkb::write(&b),
                    "{} parcel {i}, threads {threads}",
                    c["name"]
                );
            }
        }
    }
}

#[test]
fn real_utm_fixtures() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/real");
    for file in std::fs::read_dir(root).unwrap() {
        let case: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(file.unwrap().path()).unwrap()).unwrap();
        let input = geometries(&case["input"]);
        let expected_mask: Vec<bool> = case["bad"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_bool().unwrap())
            .collect();
        assert_eq!(
            validate::invalid_mask(&input),
            expected_mask,
            "{}",
            case["tile"]
        );
        for n in [1, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap();
            let mut result = input.clone();
            pool.install(|| coverage_simplify(&mut result, 5.0, Boundary::Simplify));
            for (actual, expected) in result.iter().zip(geometries(&case["output"])) {
                assert_eq!(
                    wkb::write(actual),
                    wkb::write(&expected),
                    "{}",
                    case["tile"]
                );
            }
        }
    }
}

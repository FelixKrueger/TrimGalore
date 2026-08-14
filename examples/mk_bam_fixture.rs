//! Emit the uBAM fixtures that SAM text cannot express (#415).
//!
//! BAM stores QNAME and a `Z` tag value as NUL-terminated byte fields, so either can hold a
//! newline; SAM text is line-delimited and cannot, so samtools refuses to build these. noodles'
//! record encoder rejects any read name outside `[!-?A-~]{1,254}`, so the records are packed
//! here by hand and only the BGZF framing comes from `bgzf::Writer`.
//!
//! ```sh
//! cargo run --quiet --example mk_bam_fixture -- lf_qname > test_files/ubam_lf_qname.bam
//! ```

use std::io::Write;

use noodles::bgzf;

const SEQ: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const QUAL_PHRED: u8 = 40;
const FLAG_UNMAPPED: u16 = 4;

/// 4-bit sequence alphabet, indexed by the packed nibble value.
const SEQ_CODES: &[u8] = b"=ACMGRSVTWYHKDBN";

const CASES: &[&str] = &[
    "lf_qname",
    "lf_tagvalue",
    "lf_atag",
    "atag_ok",
    "ws_qname_bulk",
];

/// Two bases per byte, high nibble first.
fn pack_seq(seq: &str) -> Vec<u8> {
    let code = |b: u8| {
        SEQ_CODES
            .iter()
            .position(|&c| c == b)
            .unwrap_or_else(|| panic!("{} is not a IUPAC base", b as char)) as u8
    };
    seq.as_bytes()
        .chunks(2)
        .map(|pair| (code(pair[0]) << 4) | pair.get(1).map_or(0, |&b| code(b)))
        .collect()
}

/// One aux field: two-char tag, `Z`, NUL-terminated value.
fn z_tag(tag: &str, value: &[u8]) -> Vec<u8> {
    let mut out = tag.as_bytes().to_vec();
    out.push(b'Z');
    out.extend_from_slice(value);
    out.push(0);
    out
}

/// One aux field: two-char tag, `A`, a single unterminated byte.
fn a_tag(tag: &str, byte: u8) -> Vec<u8> {
    let mut out = tag.as_bytes().to_vec();
    out.push(b'A');
    out.push(byte);
    out
}

/// One unmapped BAM record: `block_size` then the 32-byte fixed block, name, seq, qual, aux.
fn record(name: &[u8], aux: &[u8]) -> Vec<u8> {
    let mut name_z = name.to_vec();
    name_z.push(0);

    let mut body = Vec::new();
    body.extend_from_slice(&(-1i32).to_le_bytes()); // ref_id
    body.extend_from_slice(&(-1i32).to_le_bytes()); // pos
    body.push(u8::try_from(name_z.len()).expect("read name fits in l_read_name"));
    body.push(0); // mapq
    body.extend_from_slice(&4680u16.to_le_bytes()); // bin
    body.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op
    body.extend_from_slice(&FLAG_UNMAPPED.to_le_bytes());
    body.extend_from_slice(&(SEQ.len() as i32).to_le_bytes()); // l_seq
    body.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref_id
    body.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
    body.extend_from_slice(&0i32.to_le_bytes()); // tlen
    body.extend_from_slice(&name_z);
    body.extend_from_slice(&pack_seq(SEQ));
    body.extend_from_slice(&vec![QUAL_PHRED; SEQ.len()]);
    body.extend_from_slice(aux);

    let mut out = (body.len() as i32).to_le_bytes().to_vec();
    out.extend_from_slice(&body);
    out
}

/// The offending record sits second, so record 1 clears the sanity-check peek and the
/// failure lands in the trimming loop.
fn case_records(case: &str) -> Option<Vec<(Vec<u8>, Vec<u8>)>> {
    let clean_cb = || z_tag("CB", b"AAACCC");
    Some(match case {
        "lf_qname" => vec![
            (b"readA".to_vec(), clean_cb()),
            (b"readB\nEVIL".to_vec(), clean_cb()),
            (b"readC".to_vec(), clean_cb()),
        ],
        "lf_tagvalue" => vec![
            (b"readA".to_vec(), clean_cb()),
            (b"readB".to_vec(), z_tag("CB", b"AAA\nCCC")),
            (b"readC".to_vec(), clean_cb()),
        ],
        "lf_atag" => vec![
            (b"readA".to_vec(), a_tag("XA", b'+')),
            (b"readB".to_vec(), a_tag("XA", b'\n')),
            (b"readC".to_vec(), a_tag("XA", b'+')),
        ],
        // The acceptance twin for the `A` arm: nothing here may be refused.
        "atag_ok" => vec![
            (b"readA".to_vec(), a_tag("XA", b'+')),
            (b"readB".to_vec(), a_tag("XA", b'-')),
            (b"readC".to_vec(), a_tag("XA", b'+')),
        ],
        // #428 — enough clean records ahead of the offender that the residue is
        // large and plausible rather than one record.
        "ws_qname_bulk" => {
            let mut recs: Vec<(Vec<u8>, Vec<u8>)> = (0..5000)
                .map(|i| (format!("read{i:06}").into_bytes(), clean_cb()))
                .collect();
            recs.push((b"readEVIL with space".to_vec(), clean_cb()));
            recs.push((b"readTAIL".to_vec(), clean_cb()));
            recs
        }
        _ => return None,
    })
}

fn main() -> std::io::Result<()> {
    let case = std::env::args().nth(1).unwrap_or_default();
    let Some(records) = case_records(&case) else {
        eprintln!("usage: mk_bam_fixture <{}>", CASES.join("|"));
        std::process::exit(2);
    };

    let header_text = b"@HD\tVN:1.6\n";
    let mut writer = bgzf::Writer::new(std::io::stdout().lock());
    writer.write_all(b"BAM\x01")?;
    writer.write_all(&(header_text.len() as i32).to_le_bytes())?;
    writer.write_all(header_text)?;
    writer.write_all(&0i32.to_le_bytes())?; // n_ref
    for (name, aux) in &records {
        writer.write_all(&record(name, aux))?;
    }
    let _stdout = writer.finish()?;
    Ok(())
}

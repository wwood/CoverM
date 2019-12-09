use std::io::prelude::*;

fn main() {
    // When minimap2 indices are too big, --split-prefix is required so that @SQ
    // lines are output for samtools sort. But if --split-prefix is passed when
    // the index is not multi-part (i.e. reference has less than the number of
    // bases supplied to minimap2's -I option), then @SQ lines are output twice,
    // causing samtools to croak due to multiple references having the same
    // name.

    // The solution implemented here is to filter the minimap2 output so that
    // @SQ lines before the @PG line are removed. Other lines are printed to
    // STDOUT.

    // See https://github.com/lh3/minimap2/issues/527

    let mut is_before_pg = true;

    let stdin = std::io::stdin();
    for line_res in stdin.lock().lines() {
        match line_res {
            Ok(line) => {
                if is_before_pg {
                    if &line[0..3] == "@PG" {
                        is_before_pg = false;
                        println!("{}", line);
                    } else {
                        if &line[0..3] != "@SQ" {
                            println!("{}", line);
                        }
                    }
                } else {
                    println!("{}", line);
                }
            },
            Err(e) => {
                panic!("{}",e);
            }
        }
    }
}

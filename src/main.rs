/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

use std::process::ExitCode;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Mapped reads file
    #[arg(short, long)]
    reads: String,

    /// Genome file
    #[arg(short, long)]
    genome: String,

    /// Output file
    #[arg(short, long)]
    out: String,

    /// Be verbose
    #[arg(short, long)]
    verbose: bool,
}


fn main() -> ExitCode {

    let args = Args::parse();
    if args.verbose {
        eprintln!("[reads file={}]", args.reads);
        eprintln!("[genome file={}]", args.genome);
        eprintln!("[output file={}]", args.out);
    }

    let (names, mut chroms) = mcounts::read_fasta(&args.genome).unwrap();

    // ADS: probably should move this stuff into lib.rs?

    // make all chroms uppercase letters
    chroms.iter_mut()
        .for_each(|i| (*i).iter_mut()
                  .for_each(|j| (*j).make_ascii_uppercase()));

    mcounts::process_reads(args.verbose, args.reads,
                          args.out, names, chroms);

    ExitCode::SUCCESS
}

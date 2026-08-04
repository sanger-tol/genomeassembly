#!/usr/bin/env python
#
#    Copyright (C) 2025 Genome Research Ltd.
#
#    Author: Jim Downie <jd42@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import gzip

import click
import pyfastx


##
## fastx_format_check modified from https://github.com/lmdu/pyfastx/blob/master/pyfastxcli.py
## Author: Lianming Du
## License: MIT
##
def fastx_format_check(fastx):
    """
    Check which file format the input fastx is by checking whether it starts
    with a > or a @.
    """
    if pyfastx.gzip_check(fastx):
        fp = gzip.open(fastx, "rt")
    else:
        fp = open(fastx)

    first_char = None
    for line in fp:
        stripped = line.strip()
        if stripped:
            first_char = stripped[0]
            break

    fp.close()

    if first_char is None:
        raise RuntimeError("Error: Input file is empty!")
    elif first_char == ">":
        return "fasta"
    elif first_char == "@":
        return "fastq"
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")


@click.group()
@click.version_option(version="1.0.0", message="%(version)s")
def cli():
    """Main entry point for the tool."""
    pass


@click.command("index")
@click.argument("fastx", type=click.Path(exists=True))
def index(fastx):
    """
    Index the input fastx and print the number of sequences to stdout.
    """
    fastx_type = fastx_format_check(fastx)

    if fastx_type == "fasta":
        seq = pyfastx.Fasta(fastx, full_index=False)
    elif fastx_type == "fastq":
        seq = pyfastx.Fastq(fastx, full_index=False)
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")

    click.echo(len(seq), nl=False)


@click.command("slice")
@click.argument("fastx", type=click.Path(exists=True))
@click.argument("start", type=int)
@click.argument("end", type=int)
def slice(fastx, start, end):
    """
    Read an input fastx file and print out the slice of reads from start to end.
    """
    fastx_type = fastx_format_check(fastx)

    if fastx_type == "fasta":
        seq = pyfastx.Fasta(fastx, full_index=False)
    elif fastx_type == "fastq":
        seq = pyfastx.Fastq(fastx, full_index=False)
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")

    used_end = end
    if end > len(seq):
        click.echo(f"Warning: End position exceeds sequence length. Will only output reads until {len(seq)}", err=True)
        used_end = len(seq)

    for i in range(start, used_end):
        click.echo(seq[i].raw, nl=False)


cli.add_command(index)
cli.add_command(slice)

if __name__ == "__main__":
    cli()

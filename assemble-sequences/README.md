I found a fun application for sequence assembly.

https://pythonforbiologists.com/a-six-line-genome-assembler.html

It was implemented in Python and I decided to try to convert it into Rust to see how it performs.

So this is how we can use it:

rustc toy_assembler.rs

./toy_assembler file.fna

Note that this sequence assembler is just for fun so it won't produce the exact real consensus sequence.

# ALBATROSS implementation

## publicly AttestabLe BATched Randomness based On Secret Sharing

## Ignacio Cascudo and Bernardo David

Here, you can find a part of the implementation of the scheme ALBATROSS https://eprint.iacr.org/2020/644 over two different groups: Cyclic group and Elliptic curve group. We output the times.

#### Cyclic group
We work over a cyclic group of order q. This programm contains a mode to run the ppvss scheme of albatross, a mode to run the extraction of randomness FFTE, and a mode to run a comparison between FFTE and the extraction of randomness with a matrix filled with only 0-1.

To compile the Cyclic group code you will need the following libraries:
- GMP: https://gmplib.org/
- NTL: https://www.shoup.net/ntl/
- CryptoPP: https://www.cryptopp.com/

```
Usage:  ./albatross -p[size_of_q] [-n number_of_participants] [-h]
        ./albatross -c [-h]
        ./albatross -f[size_of_q] [-n number_of_participants] [-h]
        -p, --ppvss                   Run the ppvss scheme
                                      If the size of q isn't specified, set q of size 1024 bits
                                      Compatible with option -n
        -f, --ffte                    Run the ffte
                                      If the size of q isn't specified, set q of size 1024 bits
                                      Compatible with option -n
        -c, --comparison              Compare the two methods of extraction
                                      With n = 2048 and q of size 128
        -n, --number_of_participants  Set the number of participants n to the argument
                                      Without this option n = 1024
        -h, --help                    Display the help
Do not put space between an option and its argument, for example: ./albatross -p128 -n512
```
#### Elliptic curve  group
We work over the elliptic curve tweedledum https://eprint.iacr.org/2019/1021. The code is based on the relic-toolkit (https://github.com/relic-toolkit/relic) where this curve is implemented. In this version, only the ppvss (without the NIZK proofs) of ALBATROSS is implemented.
```
Usage:  ./albatross [n]
        if n isn't specified, n = 2048
```

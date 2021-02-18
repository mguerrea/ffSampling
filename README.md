# ffSampling

## TEST mode

`make TEST=true && ./ffsampling`

## Standard mode

`make`

(It may be required to clean object files before with `make clean` if the project was compiled in TEST mode)

`Usage: ./ffsampling [dim] [-s|-v] [message]`

-s: performs signature of message

-v: verifies signature of message

Private keys are hardcoded because the only purpose of this program is to provide a demonstration for the ffSampling algorithm.

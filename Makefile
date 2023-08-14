# Specify the target executable and the source files needed to build it
my_app: symnmf.o symnmfmodule.o symnmf.h
    gcc -o my_app symnmf.o symnmfmodule.o
# Specify the object files that are generated from the corresponding source files
symnmf.o: symnmf.c
    gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.c -lm

symnmfmodule.o: symnmfmodule.c
    gcc -c foo.c

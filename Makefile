NAME = ffsampling

SRC_DIR := srcs
OBJ_DIR := ./obj
INC_DIR := includes

LIBS = -lgmp -lm
CC = clang

C_FILES = samplerz.c tables.c random.c roots.c fft.c utils.c pol_op.c \
	keygen.c params.c ffsampling.c ntt.c

TEST ?= false

ifeq ($(TEST), false)
	C_FILES += main.c
else
	C_FILES += test.c KAT.c
	C_FLAGS = -DTEST
endif

H_FILES = ffsampling.h

### MAIN AND SUB FILES ###
O_FILES = $(C_FILES:.c=.o)

### Full Paths ###
SRC = $(addprefix $(SRC_DIR)/,$(C_FILES))
OBJ = $(addprefix $(OBJ_DIR)/,$(O_FILES))
INC = $(addprefix $(INC_DIR)/,$(H_FILES))

INC_FLAGS = $(addprefix -I ,$(INC_DIR))

all: $(NAME)

$(NAME): $(OBJ_DIR) $(OBJ)
	@$(CC) -o $@ $(OBJ) $(LIBS)
	@echo "[$(NAME)] Compilation success"

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INC)
	@$(CC) -c $< -o $@ $(INC_FLAGS) $(C_FLAGS)
	@echo "[$(NAME)] $@ created"

clean:
	@/bin/rm -rf $(OBJ)
	@/bin/rm -rf $(OBJ_DIR)
	@echo "[$(NAME)] .o files deleted$()$(RESET)"

fclean: clean
	@/bin/rm -f $(NAME)
	@/bin/rm -f $(LINKNAME)
	@echo  "[$(NAME)] executable file deleted$(RESET)"

re: fclean $(NAME)

.PHONY: all, clean, fclean, re
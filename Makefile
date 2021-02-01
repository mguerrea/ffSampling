NAME = ffsampling

SRC_DIR := .
OBJ_DIR := ./obj
INC_DIR := .

LIBS = -lgmp -lm
CC = clang

C_FILES = samplerz.c tables.c random.c roots.c fft.c utils.c fft_op.c \
	keygen.c params.c

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
	@echo "\n$(BLU)[$(NAME)]$(GRN) Compilation success$(RESET)"

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INC)
	$(CC) -c $< -o $@ $(INC_FLAGS) $(C_FLAGS)
	@echo "\r$(ERASE)$(BLU)[$(NAME)]$(RESET) $@ created\c"

clean:
	@/bin/rm -rf $(OBJ)
	@/bin/rm -rf $(OBJ_DIR)
	@echo "$(BLU)[$(NAME)]$(RED) .o files deleted$()$(RESET)"

fclean: clean
	@/bin/rm -f $(NAME)
	@/bin/rm -f $(LINKNAME)
	@echo  "$(BLU)[$(NAME)]$(RED) executable file deleted$(RESET)"

re: fclean $(NAME)

.PHONY: all, clean, fclean, re
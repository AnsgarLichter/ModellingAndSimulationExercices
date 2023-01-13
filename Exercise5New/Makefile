TARGET := dynamischessystem
HEADER := $(TARGET).h
SOURCES := $(TARGET).c
OBJS := $(SOURCES:.c=.o)
CFLAGS := -Wall -O1
LIBS := -lm
CC := $(shell which gcc)

all: $(TARGET)

$(TARGET): $(OBJS) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean :
	rm $(TARGET) $(TARGET).o
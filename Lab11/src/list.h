#ifndef LIST_H 
#define LIST_H

#include <stdio.h>
#include <stdlib.h>

typedef struct node node;
struct node {
    char data;
    node *next;
};

void generate_random_list(node **head, int n);
void print_list(node *head);
void print_list_reverse(node *head);

int search_list(node *head, char target);
char get_element(node *head, int index);

int length(node *head);
char pop(node **head);
void push(node **head, char value);
void delete_list(node **head);

#endif

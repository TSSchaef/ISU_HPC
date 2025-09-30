#include "list.h"

int search_list(node *head, char target){
    int i;
    for(i = 0; head != NULL; head = head->next, i++){
        if(head->data == target) return i;
    }
    return -1;
}

char get_element(node *head, int index){
    int i;
    for(i = 0; head != NULL; head = head->next, i++){
        if(i == index) return head->data;
    }
    return '\0';
}

int length(node *head){
    int i;
    for(i = 0; head != NULL; head = head->next, i++);
    return i;
}

char pop(node **head){
    if(*head == NULL) return '\0';

    node *temp = *head;
    char value = temp->data;
    *head = (*head)->next;
    free(temp);
    return value;
}

void push(node **head, char value){
    node *new_node = (node *)malloc(sizeof(node));
    new_node->data = value;
    new_node->next = *head;
    *head = new_node;
}

void delete_list(node **head){
    while(*head != NULL){
        pop(head);
    }
}

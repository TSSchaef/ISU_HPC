#include "list.h"

void generate_random_list(node **head, int n){
    int i;
    for(i = 0; i < n; i++){
        char value;
        // Generate random uppercase or lowercase letter
        if( rand() % 2 == 0 )
            //ASCII arithmetic
            value = 'A' + rand() % 26;
        else {
            value = 'a' + rand() % 26;
        }

        push(head, value);
    }
}

void print_list(node *head){
    for(; head != NULL; head = head->next){
        printf("%c -> ", head->data);
    }
    printf("NULL\n");
}

// Print list in reverse using right recursion
void print_list_reverse_recursive(node *head){
    if(head == NULL) {
        printf("NULL");
        return;
    }
    print_list_reverse_recursive(head->next);
    printf(" <- %c", head->data);
}

void print_list_reverse(node *head){
    print_list_reverse_recursive(head);
    printf("\n");
}



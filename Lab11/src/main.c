#include <stdio.h>
#include "list.h"


node *head = NULL;


int main(int argc, char *argv[]) {
    printf("Generating random list:\n");
    generate_random_list(&head, 10);
    print_list(head);

    printf("Printing list in reverse order:\n");
    print_list_reverse(head);

    printf("Pushing 's' to the front of the list:\n");
    push(&head, 's');
    print_list(head);

    printf("Popping the front element of the list: %c \n", pop(&head));
    printf("Popping the front element of the list: %c \n", pop(&head));

    printf("Current list:\n");
    print_list(head);

    printf("Searching for 'a' in the list: %d (-1 means not in the list)\n", 
           search_list(head, 'a'));
    printf("Searching for 'K' in the list: %d (-1 means not in the list)\n", 
           search_list(head, 'K'));

    printf("Getting element at index 3: %c \n", get_element(head, 3));

    printf("The length of the list is: %d\n", length(head));

    printf("Deleteing the entire list.\n");
    delete_list(&head);
    printf("Remaining list (should be empty):\n");
    print_list(head);

    return 0;
}

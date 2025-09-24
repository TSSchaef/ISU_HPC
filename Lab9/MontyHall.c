#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    int carDoor, userChoice, montyOpens, finalChoice;
    char switchChoice;
    char playAgain = 'y';

    srand(time(NULL));
    carDoor = rand() % 3 + 1;

    printf("Welcome to the Monty Hall game!\n");
    printf("There are 3 doors. Behind one is a car, behind the others are goats.\n");


    printf("Pick a door (1, 2, or 3): ");
    scanf("%d", &userChoice);

    // Monty opens a door with a goat that is not the user's choice or the car
    do {
        montyOpens = rand() % 3 + 1;
    } while (montyOpens == userChoice || montyOpens == carDoor);

    printf("Monty opens door %d to reveal a goat!\n", montyOpens);

    printf("Do you want to switch doors? (y/n): ");
    scanf(" %c", &switchChoice);

    if (switchChoice == 'y' || switchChoice == 'Y') {
        // Switch choice to the only remaining unopened door
        for (int i = 1; i <= 3; i++) {
            if (i != userChoice && i != montyOpens) {
                finalChoice = i;
                break;
            }
        }
    } else {
        finalChoice = userChoice;
    }

    if (finalChoice == carDoor) {
        printf("Congratulations! You won the car behind door %d!\n", finalChoice);
    } else {
        printf("Sorry, you got a goat behind door %d.\n", finalChoice);
    }

    printf("The car was behind door %d.\n", carDoor);

    printf("Input q to quit, enter any other character to play again!");
    scanf(" %c", &playAgain);

    if (playAgain != 'q' && playAgain != 'Q') {
        main(); // Restart the game
    }

    return 0;
}

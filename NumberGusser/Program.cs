using System;

namespace NumberGuesser
{
    class Program
    {
        private static void Main(string[] args)
        {
            // App Info
            GetAppInfo();
            GreetUser();

            // Initialize random number generator
            Random random = new Random();

            while (true)
            {
                // Generate random number
                int correctNumber = random.Next(1, 10);
                int guess = 0;

                Console.WriteLine("Guess a number between 1 and 10");

                // Loop until the guess is correct
                while (guess != correctNumber)
                {
                    string input = Console.ReadLine();

                    // Validate input
                    if (!int.TryParse(input, out guess))
                    {
                        PrintColorMessage(ConsoleColor.Red, "Please enter an actual number");
                        continue;
                    }

                    // Compare guess with the correct number
                    if (guess != correctNumber)
                    {
                        PrintColorMessage(ConsoleColor.Red, "Wrong Number, please try again");
                    }
                }

                // Success message
                PrintColorMessage(ConsoleColor.Yellow, "You are correct!");

                // Ask if the player wants to play again
                Console.WriteLine("Play again? [Y or N]");
                string answer = Console.ReadLine().ToUpper();

                if (answer == "Y")
                {
                    continue;
                }
                else if (answer == "N")
                {
                    break;
                }
                else
                {
                    break;
                }
            }
        }

        // Method to display application info
        private static void GetAppInfo()
        {
            string appName = "Number Guesser";
            string appVersion = "1.0.0";
            string appAuthor = "Stavros Zymvragoudakis";

            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine("{0}: Version {1} by {2}", appName, appVersion, appAuthor);
            Console.ResetColor();
        }

        // Method to greet the user
        private static void GreetUser()
        {
            Console.WriteLine("What is your name?");
            string inputName = Console.ReadLine();
            Console.WriteLine("Hello {0}, let's play a game



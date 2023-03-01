import random

# set the number of persons, meetings and runs
number_of_persons = 4
number_of_meetings = 4
number_of_runs = 1000
contact_probability = 1/2

# dictionary to store the number of contacts per meeting in regard to atendees
contacts_per_meeting = {0: 0, 1: 0, 2: 2, 3: 6, 4: 12}


sum_of_contacts = 0
# three for loops to iterate
for i in range(1, number_of_runs+1):
    for j in range(1, number_of_meetings+1):
        number_of_persons_in_meeting = 0
        for k in range(1, number_of_persons+1):
            # if the random number is less than 0.5, the person will be at the meeting
            if random.random() < contact_probability:
                number_of_persons_in_meeting = number_of_persons_in_meeting+1
        # add the number of contacts to the sum
        sum_of_contacts = sum_of_contacts + contacts_per_meeting[number_of_persons_in_meeting]

# calculate the average number of contacts
average_number_of_contacts_per_meeting = sum_of_contacts / (number_of_runs * number_of_meetings)
average_number_of_contacts_per_run = sum_of_contacts / number_of_runs

print("The average number of contacts per meeting is: ", average_number_of_contacts_per_meeting)
print("The average number of contacts per run is: ", average_number_of_contacts_per_run)

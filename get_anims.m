function anims = get_anims(abl_type)
% gets all relevant animals so you can harmonize across analyses as you add animals

    switch abl_type
        case 'touch'
            anims = {'250220','257220','258836','271211', '281915','274424','272761', '274577','275801'};    

        case 'whisking'
            anims = {'278937','278939', '281915', '275798', '278288', '276013', '278759'};

        case 'silent'
            anims = {'257218', '258836', '271211', '278937', '278939', '275801', '275798', '276013'}; 
    end

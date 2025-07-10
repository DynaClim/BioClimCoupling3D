#!/bin/bash

# Optionnel : intervalle pour attendre que les jobs "will never run" apparaissent
sleep 10

# Récupérer tous les jobs de l'utilisateur dans l'état "DEPENDENCY"
jobs_to_cancel=$(sacct --format=JobID,State --noheader | awk '$2=="DEPENDENCY" {print $1}')

if [ -z "$jobs_to_cancel" ]; then
    echo "Aucun job à annuler."
else
    echo "Annulation des jobs suivants :"
    echo "$jobs_to_cancel"
    scancel $jobs_to_cancel
fi

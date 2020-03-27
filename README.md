# epidemilogy
Compartmental models may be used to predict properties of how a disease spreads, for example the prevalence (total number of infected) or the duration of an epidemic. Also, the model allows for understanding how different situations may affect the outcome of the epidemic, e.g., what the most efficient technique is for issuing a limited number of vaccines in a given population; 
see https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model for more details.
This is based on the really nice explanations in https://towardsdatascience.com/social-distancing-to-slow-the-coronavirus-768292f04296 and many code comments are taken from there.

# Running the app

First a small webserver has to be started. To do so, type

```bash
python -m http.server
```

Then open your browser on http://localhost:8000
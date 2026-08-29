These are two expert-level plant biology multiple-choice datasets with identical schemas, used to evaluate LLM performance before and after integrating the PlantGeneCopilot MCP. All questions are is_expert=True (expert-annotated), with 3 options each.
Shared fields (per question)
question / area (fine-grained category) / plant_species (species list) / options (3-choice) / source + doi + source_journal (provenance) / Year / Citations / answer (index of the correct option) / normalized_plant_species / normalized_area / is_expert
1. filtered_species.json — Target plant (308 questions)
- Focuses on model organisms, covering Arabidopsis thaliana (dominant), tomato, rice, maize, Marchantia, etc.; 12 questions involve multiple species.
- Normalized species distribution: Model Organisms 285, Solanaceae & Relatives 17, Cereal Grains 4, Woody Perennials & Trees 2.
- Subject-area distribution: ENVIRONMENT 77, GENE REGULATION 73, GROWTH AND DEVELOPMENT 54, GENOME AND GENOMICS 29, PHYSIOLOGY AND METABOLISM 20, CELL BIOLOGY AND CELL SIGNALING 20, HORMONES 18, PLANT BIOTECHNOLOGY 9, EVOLUTION 8.
- Years 1994–2024, mostly recent, highly-cited literature.
2. other_species.json — Non-target plant (257 questions)
- Focuses on non-model / species-specific plants, with a large share of Non-specific (cross-species/general) questions (167).
- Normalized species: Non-specific 167, Model Organisms 31, Legumes 16, Woody Perennials & Trees 15, Solanaceae & Relatives 13, Other Herbaceous Crops, Spices, Fibers & Weeds 9, Cereal Grains 6.
- Subject-area distribution: GENE REGULATION 78, ENVIRONMENT 75, PLANT BIOTECHNOLOGY 26, GROWTH AND DEVELOPMENT 23, HORMONES 18, EVOLUTION 16, GENOME AND GENOMICS 12, PHYSIOLOGY AND METABOLISM 6, CELL BIOLOGY AND CELL SIGNALING 3.
- Same year range 1994–2024.
Together, these form the target-plant vs. non-target/species-specific plant evaluation benchmark (i.e., the target / non-target accuracy data reported earlier).

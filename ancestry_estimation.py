import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix

# Load the .evec file
evec_file = 'CASES_and_CONTROLS_1kg_merge.pca.evec_RENAMED'
evec_data = pd.read_csv(evec_file, delim_whitespace=True, header=None)

# Rename columns (assuming last column is population or case/control labels)
n_pcs = evec_data.shape[1] - 2  # Number of PCs (subtract 2 for SampleID and Pop/Label columns)
evec_data.columns = ['SampleID'] + [f'PC{i}' for i in range(1, n_pcs + 1)] + ['Pop']

# Separate 1000 Genomes reference samples (those with ancestry labels)
known_samples = evec_data[~evec_data['Pop'].isin(['CASE', 'CONTROL'])]

# Separate CASE samples
case_samples = evec_data[evec_data['Pop'] == 'CASE']

# Separate CONTROL samples
control_samples = evec_data[evec_data['Pop'] == 'CONTROL']

# Features: Use the PCs (PC1, PC2, ..., PCn)
X = known_samples[[f'PC{i}' for i in range(1, n_pcs + 1)]]

# Target: Population labels (e.g., AFR, EUR, EAS, SAS, AMR)
y = known_samples['Pop']

# Split into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize the Random Forest classifier
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
rf_model.fit(X_train, y_train)

# Evaluate the model on the test set
y_pred = rf_model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

print(accuracy)
print(conf_matrix)

# Extract PCs for the CASE samples
X_case = case_samples[[f'PC{i}' for i in range(1, n_pcs + 1)]]

# Predict ancestry for CASE samples
predicted_case_ancestry = rf_model.predict(X_case)

# Get prediction probabilities for CASE samples (for certainty)
predicted_case_proba = rf_model.predict_proba(X_case)

# Calculate certainty and add it to the CASE samples dataframe
case_samples['Certainty'] = [
    predicted_case_proba[i][rf_model.classes_.tolist().index(predicted_case_ancestry[i])] 
    for i in range(len(predicted_case_proba))
]

# Add predicted ancestry to the CASE samples dataframe
case_samples['Predicted Ancestry'] = predicted_case_ancestry

# Output CASE samples with predicted ancestry and certainty
case_samples.to_csv('cases_estimated_ancestry.csv', index=False)

# Extract PCs for the CONTROL samples
X_control = control_samples[[f'PC{i}' for i in range(1, n_pcs + 1)]]

# Predict ancestry for CONTROL samples
predicted_control_ancestry = rf_model.predict(X_control)

# Get prediction probabilities for CONTROL samples (for certainty)
predicted_control_proba = rf_model.predict_proba(X_control)

# Calculate certainty and add it to the CONTROL samples dataframe
control_samples['Certainty'] = [
    predicted_control_proba[i][rf_model.classes_.tolist().index(predicted_control_ancestry[i])] 
    for i in range(len(predicted_control_proba))
]

# Add predicted ancestry to the CONTROL samples dataframe
control_samples['Predicted Ancestry'] = predicted_control_ancestry

# Combine case_samples and control_samples DataFrames
combined_samples = pd.concat([case_samples, control_samples])

# Write combined DataFrame to a CSV file with certainty
combined_samples.to_csv('cases_and_controls_estimated_ancestry.csv', index=False)
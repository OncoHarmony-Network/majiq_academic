import numpy as np
from sklearn.model_selection import StratifiedKFold

n_samples = 4
X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
y = np.array([10, 11, 12, 13])
skf = StratifiedKFold(n_splits=2)
print(skf.get_n_splits(X, y))

print(skf)

StratifiedKFold(n_splits=2, random_state=None, shuffle=False)

for train_index, test_index in skf.split(X, y):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

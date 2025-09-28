from utils.tgtinputhandler import *
from pairgen import *
from filtering import *
from featuresforval import *
from ranking.rankingtwo import *
import csv




primer_features = FeatureCalculating("src/weight_val.csv")
print("got the dataset")
primer_features.compute_all_features()
print("calculated features")
primer_features.to_csv("val_with_features.csv")

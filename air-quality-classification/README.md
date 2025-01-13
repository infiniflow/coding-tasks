#### Task description

1. Prediction of air quality (Good/Moderate/Hazardous/Poor) based on the data of the given weather indicators (Temperature/Humidity/PM2.5/PM10/NO2/SO2/CO/Proximity to Industrial Areas/Population Density).
2. Give the model training code (Python 3.11), and the trained model.
3. Support with docker image build and `docker compose up` one click to start the backend service (Python 3.11).
4. You can use the browser to open the page and then enter the value of each indicator (allowed to be empty), submit and return to the prediction results.

#### Dataset

You need to use following [dataset](./updated_pollution_dataset.csv) to train your model.

#### Submission

You need to package your submission into `aqc.tar.gz` and send it to `kai.hu@infiniflow.ai`.

#### Evalutation

1. We are going to use following instruction to start your service:

```shell
tar -xvf aqc.tar.gz
cd aqc
docker build -t evaluation/aqc .
cd docker
docker compose up -d
```

2. The page for inputting indicators can be opened by opening a browser and going to the following address: http://host:1111/aqc.

3. Enter the indicators for each dimension and click OK to return the predicted results (Good/Moderate/Hazardous/Poor).

4. Access to the container allows the execution of training scripts, as recommended below:

   ```shell
   docker exec -it aqc /bin/bash
   python train.py dataset.csv
   #print F1 score.
   ```

   

@description('Location for all resources.')
param location string = resourceGroup().location

@description('Base name for the apps')
param appBaseNae string

@description('The image to use for zookeeper')
param zookeeperImage string = 'confluentinc/cp-zookeeper:7.6.0'

@description('The image to use for kafka')
param kafkaImage string = 'confluentinc/cp-kafka:7.6.0'

@description('The image to use for kafka-connect')
param kafkaConnectImage string = 'confluentinc/cp-kafka-connect:7.6.0'

@description('Resource ID of the managed environment for Azure Container Apps.')
param containerAppEnvId string

// Zookeeper settings
var zookeeperClientPort string = '2181'
var zookeeperTickTime string = '2000'

// Kafka settings
var kafkaMainPort string = '9092'
var kafkaAltPort string = '29092'


// Kafka Connect settings


// Kafka stack container app
resource containerApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: '${appBaseNae}-kafka-stack'
  location: location
  properties: {
    managedEnvironmentId: containerAppEnvId
    template: {
      containers: [
        // Zookeeper
        {
          name: '${appBaseNae}-zookeeper'
          image: zookeeperImage
          env: [
            {
              name: 'ZOOKEEPER_CLIENT_PORT'
              value: zookeeperClientPort
            }
            {
              name: 'ZOOKEEPER_TICK_TIME'
              value: zookeeperTickTime
            }
          ]
        }
        // Kafka
        {
          name: '${appBaseNae}-kafka'
          image: kafkaImage
          env: [
            {
              name: 'KAFKA_BROKER_ID'
              value: '1'
            }
            {
              name: 'KAFKA_ZOOKEEPER_CONNECT'
              value: 'localhost:${zookeeperClientPort}'
            }
            {
              name: 'KAFKA_LISTENER_SECURITY_PROTOCOL_MAP'
              value: 'PLAINTEXT:PLAINTEXT,PLAINTEXT_HOST:PLAINTEXT'
            }
            {
              name: 'KAFKA_ADVERTISED_LISTENERS'
              value: 'PLAINTEXT://localhost:${kafkaAltPort},PLAINTEXT_HOST://localhost:${kafkaMainPort}'
            }
            {
              name: 'KAFKA_LISTENERS'
              value: 'PLAINTEXT://0.0.0.0:${kafkaAltPort},PLAINTEXT_HOST://0.0.0.0:${kafkaMainPort}'
            }
            {
              name: 'KAFKA_INTER_BROKER_LISTENER_NAME'
              value: 'PLAINTEXT'
            }
            {
              name: 'KAFKA_OFFSETS_TOPIC_REPLICATION_FACTOR'
              value: '1'
            }
          ]
        }
        // Kafka-connect
        {
          name: '${appBaseNae}-kafka-connect'
          image: kafkaConnectImage
          env: [
            {
              name: 'CONNECT_BOOTSTRAP_SERVERS'
              value: 'localhost:${kafkaAltPort}'
            }
            {
              name: 'CONNECT_REST_PORT'
              value: '8083'
            }
            {
              name: 'CONNECT_GROUP_ID'
              value: 'kafka-connect-group'
            }
            {
              name: 'CONNECT_CONFIG_STORAGE_TOPIC'
              value: 'connect-configs'
            }
            {
              name: 'CONNECT_OFFSET_STORAGE_TOPIC'
              value: 'connect-offsets'
            }
            {
              name: 'CONNECT_STATUS_STORAGE_TOPIC'
              value: 'connect-status'
            }
            {
              name: 'CONNECT_KEY_CONVERTER'
              value: 'org.apache.kafka.connect.storage.StringConverter'
            }
            {
              name: 'CONNECT_VALUE_CONVERTER'
              value: 'org.apache.kafka.connect.json.JsonConverter'
            }
            {
              name: 'CONNECT_VALUE_CONVERTER_SCHEMAS_ENABLE'
              value: 'false'
            }
            {
              name: 'CONNECT_KEY_CONVERTER_SCHEMAS_ENABLE'
              value: 'false'
            }
            {
              name: 'CONNECT_INTERNAL_KEY_CONVERTER'
              value: 'org.apache.kafka.connect.json.JsonConverter'
            }
            {
              name: 'CONNECT_INTERNAL_VALUE_CONVERTER'
              value: 'org.apache.kafka.connect.json.JsonConverter'
            }
            {
              name: 'CONNECT_LOG4J_LOGGERS'
              value: 'org.reflections=ERROR'
            }
            {
              name: 'CONNECT_PLUGIN_PATH'
              value: '/usr/share/java,/etc/kafka-connect/jars'
            }
            {
              name: 'CONNECT_REST_ADVERTISED_HOST_NAME'
              value: 'kafka-connect'
            }
          ]
        }
      ]
    }
  }
}

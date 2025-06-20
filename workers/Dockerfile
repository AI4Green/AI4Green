# Use official Python image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install Poetry
RUN pip install poetry==1.8.3

# Copy project files into the container
COPY . .

# Install dependencies via Poetry
RUN poetry install --no-dev

# Command to run the worker script
CMD ["poetry", "run", "python", "kafka_worker.py"]

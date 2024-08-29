# 29/08/24  -   This docker file sets up a light-weight version of AI4Green in a container, but for some reasone wont
#               connect to a postrgres database configured in a seperate container. We've removed the option to install AI4Green
#               via docker from the REDME.md for now

# Use the official slim Python image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    libxrender1 \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
RUN pip install poetry==1.8.3

# Set the working directory
WORKDIR /app

# Copy only the pyproject.toml and poetry.lock files first
COPY pyproject.toml poetry.lock ./

# Install Python dependencies using Poetry
RUN poetry config virtualenvs.create false \
    && poetry install --no-root --no-interaction --no-ansi

# Copy the entire project to the working directory
COPY . .

# Expose the necessary port
EXPOSE 80

# Run the application
CMD ["python", "Webapp/app.py"]

params.ncbi_api_key = '<Your API key here>'

Channel
  .fromSRA(['SRP073307'], apiKey: params.ncbi_api_key)
  .view()
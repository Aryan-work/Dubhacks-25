import { StatsigClient } from '@statsig/js-client';
import { StatsigSessionReplayPlugin } from '@statsig/session-replay';
import { StatsigAutoCapturePlugin } from '@statsig/web-analytics';

class StatsigService {
  constructor() {
    this.client = null;
    this.isInitialized = false;
    this.clientKey = "client-WlGqgDoQ3OgrjFxQPFaQhOi1LebsJ3kKby9NgdYUCgX";
  }

  async initialize(userID = null) {
    if (this.isInitialized) {
      return;
    }

    try {
      // Generate a unique user ID if none provided
      const userId = userID || this.generateUserId();
      
      this.client = new StatsigClient(
        this.clientKey,
        { userID: userId },
        {
          plugins: [
            new StatsigSessionReplayPlugin(),
            new StatsigAutoCapturePlugin(),
          ],
        }
      );

      await this.client.initializeAsync();
      this.isInitialized = true;
      
      console.log('✅ Statsig initialized successfully with session replay and auto-capture');
      
      // Log initialization event
      this.logEvent('statsig_initialized', {
        user_id: userId,
        timestamp: new Date().toISOString(),
        features: ['session_replay', 'auto_capture', 'web_analytics']
      });
      
    } catch (error) {
      console.error('❌ Failed to initialize Statsig:', error);
    }
  }

  generateUserId() {
    // Generate a unique user ID for anonymous users
    return `user_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
  }

  async checkGate(gateName, userID = null) {
    if (!this.isInitialized) {
      await this.initialize(userID);
    }
    
    try {
      return await this.client.checkGate(gateName);
    } catch (error) {
      console.error(`Error checking gate ${gateName}:`, error);
      return false;
    }
  }

  async getConfig(configName, userID = null) {
    if (!this.isInitialized) {
      await this.initialize(userID);
    }
    
    try {
      return await this.client.getConfig(configName);
    } catch (error) {
      console.error(`Error getting config ${configName}:`, error);
      return {};
    }
  }

  async getExperiment(experimentName, userID = null) {
    if (!this.isInitialized) {
      await this.initialize(userID);
    }
    
    try {
      return await this.client.getExperiment(experimentName);
    } catch (error) {
      console.error(`Error getting experiment ${experimentName}:`, error);
      return {};
    }
  }

  logEvent(eventName, value = null, metadata = {}) {
    if (!this.isInitialized) {
      console.warn('Statsig not initialized, cannot log event:', eventName);
      return;
    }
    
    try {
      this.client.logEvent(eventName, value, metadata);
      console.log(`📊 Event logged: ${eventName}`, { value, metadata });
    } catch (error) {
      console.error(`Error logging event ${eventName}:`, error);
    }
  }

  // Molecular analysis specific events
  logMolecularAnalysis(molecularData, analysisType, results) {
    this.logEvent('molecular_analysis_started', {
      input_type: analysisType,
      molecular_data: molecularData.substring(0, 100), // Truncate for privacy
      timestamp: new Date().toISOString()
    });
  }

  logQuantumSimulation(simulationConfig, results) {
    this.logEvent('quantum_simulation_completed', {
      config: simulationConfig,
      results_summary: {
        homo_lumo_gap: results?.homo_lumo_gap,
        binding_energy: results?.binding_energy,
        molecular_analogs_count: results?.molecular_analogs?.length || 0
      },
      timestamp: new Date().toISOString()
    });
  }

  logUserInteraction(interaction, component, metadata = {}) {
    this.logEvent('user_interaction', {
      interaction,
      component,
      metadata,
      timestamp: new Date().toISOString()
    });
  }

  logError(error, context = {}) {
    this.logEvent('error_occurred', {
      error_message: error.message,
      error_stack: error.stack,
      context,
      timestamp: new Date().toISOString()
    });
  }

  // Feature flags for molecular analysis platform
  async shouldShowAdvancedVisualization(userID = null) {
    return await this.checkGate('advanced_visualization_gate', userID);
  }

  async shouldEnableQuantumMode(userID = null) {
    return await this.checkGate('quantum_mode_gate', userID);
  }

  async getAnalysisAlgorithm(userID = null) {
    const config = await this.getConfig('analysis_algorithm_config', userID);
    return config.get('algorithm', 'standard');
  }

  async getVisualizationSettings(userID = null) {
    const config = await this.getConfig('visualization_settings', userID);
    return {
      defaultRepresentation: config.get('default_representation', 'ballstick'),
      enableRotation: config.get('enable_rotation', true),
      showLabels: config.get('show_labels', true)
    };
  }

  // Session replay and analytics
  startSession() {
    this.logEvent('session_started', {
      timestamp: new Date().toISOString(),
      user_agent: navigator.userAgent,
      screen_resolution: `${window.screen.width}x${window.screen.height}`
    });
  }

  endSession() {
    this.logEvent('session_ended', {
      timestamp: new Date().toISOString()
    });
  }

  // Cleanup
  async shutdown() {
    if (this.client) {
      await this.client.shutdown();
      this.isInitialized = false;
    }
  }
}

// Create singleton instance
const statsigService = new StatsigService();

export default statsigService;

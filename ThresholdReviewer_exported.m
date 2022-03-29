classdef ThresholdReviewer_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        Panel                         matlab.ui.container.Panel
        GridLayout2                   matlab.ui.container.GridLayout
        GridLayout4                   matlab.ui.container.GridLayout
        SaveTruePositiveEventsButton  matlab.ui.control.Button
        GridLayout6                   matlab.ui.container.GridLayout
        GridLayout8                   matlab.ui.container.GridLayout
        SharpwaveCheckBox             matlab.ui.control.CheckBox
        TruePositiveSharpWaveLabel    matlab.ui.control.Label
        GridLayout5                   matlab.ui.container.GridLayout
        GridLayout7                   matlab.ui.container.GridLayout
        RippleCheckBox                matlab.ui.control.CheckBox
        TruePositiveRippleLabel       matlab.ui.control.Label
        GridLayout3                   matlab.ui.container.GridLayout
        EventsSlider                  matlab.ui.control.Slider
        NextButton                    matlab.ui.control.Button
        LastButton                    matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        fn = 600;
        waveforms;
        oscil_table;
        dataset_path;
        events_number;
        window_length = 3601;
        true_positive_ripples = [];
        true_positive_sharpwaves = [];
    end
    
    methods (Access = private)
        function update(app)
            event_num = app.EventsSlider.Value;
            
            if ismember(event_num, app.true_positive_ripples)
                app.RippleCheckBox.Value = 1;
            else
                app.RippleCheckBox.Value = 0;
            end
            if ismember(event_num, app.true_positive_sharpwaves)
                app.SharpwaveCheckBox.Value = 1;
            else
                app.SharpwaveCheckBox.Value = 0;
            end
            
            x = (1:app.window_length)./app.fn;
            
            belo = app.waveforms.HPCbelo;
            belo_sig = belo(:,event_num);
            belo_width = max(belo_sig) - min(belo_sig);
            
            pyra = app.waveforms.HPCpyra;
            pyra_sig = pyra(:,event_num);
            pyra_sig = pyra_sig + belo_width;
            
            title(app.UIAxes, strcat("Event ", int2str(event_num)))
            cla(app.UIAxes)
            plot(app.UIAxes, x, belo_sig, 'Color', 'black')
            hold(app.UIAxes,'on')
            plot(app.UIAxes, x, pyra_sig, 'Color', 'blue')
            
            app.UIAxes.YLim = [min(belo_sig)-500, max(pyra_sig)+500];
            app.UIAxes.XLim = [0, app.window_length/app.fn];
            hold(app.UIAxes,'on')
            txt = app.oscil_table.Type(event_num);
            text(app.UIAxes, app.window_length/app.fn-0.5, max(pyra_sig)+250, txt)

        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %D:\Dev\MATLAB\results_sleep_oscil_tools_batch2.1\dataset201.mat
            [file,path] = uigetfile();
            if file == 0
                delete(app); 
            else
                app.dataset_path = fullfile(path,file);
            
                table_name = 'grouped_oscil_table';
                oscil_table_ = load(app.dataset_path, table_name);
                app.oscil_table = oscil_table_.grouped_oscil_table;
                
                waveforms_name = 'grouped_wave_forms';
                waveforms_ = load(app.dataset_path, waveforms_name);
                app.waveforms = waveforms_.grouped_wave_forms;
                
                app.events_number = height(app.oscil_table);
                app.EventsSlider.Limits = [1 app.events_number];
                
                app.SaveTruePositiveEventsButton.Enable = 'on';
                app.SharpwaveCheckBox.Enable = 'on';
                app.RippleCheckBox.Enable = 'on';
                app.EventsSlider.Enable = 'on';
                app.NextButton.Enable = 'on';
                app.LastButton.Enable = 'on';
                
                app.update()
            end
        end

        % Button pushed function: NextButton
        function NextButtonPushed(app, event)
            if app.EventsSlider.Value < app.events_number
                app.EventsSlider.Value = app.EventsSlider.Value + 1;
                app.update()
            end
        end

        % Button pushed function: LastButton
        function LastButtonPushed(app, event)
            if app.EventsSlider.Value > 1
                app.EventsSlider.Value = app.EventsSlider.Value - 1;
                app.update()
            end
        end

        % Value changed function: EventsSlider
        function EventsSliderValueChanged(app, event)
            value = round(app.EventsSlider.Value);
            app.EventsSlider.Value = value;
            app.update()
        end

        % Value changed function: RippleCheckBox
        function RippleCheckBoxValueChanged(app, event)
            value = app.RippleCheckBox.Value;
            if value
                app.true_positive_ripples = [app.true_positive_ripples app.EventsSlider.Value];
            else
                app.true_positive_ripples(app.true_positive_ripples == app.EventsSlider.Value) = [];
            end
        end

        % Value changed function: SharpwaveCheckBox
        function SharpwaveCheckBoxValueChanged(app, event)
            value = app.SharpwaveCheckBox.Value;
            if value
                app.true_positive_sharpwaves = [app.true_positive_sharpwaves app.EventsSlider.Value];
            else
                app.true_positive_sharpwaves(app.true_positive_sharpwaves == app.EventsSlider.Value) = [];
            end
        end

        % Button pushed function: SaveTruePositiveEventsButton
        function SaveTruePositiveEventsButtonPushed(app, event)
            [fPath, fName, fExt] = fileparts(app.dataset_path);
            file = strcat(fName,'TruePositives',fExt);
            path = fPath;
            true_positives_file = fullfile(path,file);
            ripples = app.true_positive_ripples;
            sharpwaves = app.true_positive_sharpwaves;
            save(true_positives_file,'ripples', 'sharpwaves');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'2x', '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Event')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 1;

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Layout.Row = 2;
            app.Panel.Layout.Column = 1;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel);
            app.GridLayout2.ColumnWidth = {'1x'};

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout2);
            app.GridLayout3.ColumnWidth = {'1x', '8x', '1x'};
            app.GridLayout3.RowHeight = {'1x'};
            app.GridLayout3.Layout.Row = 2;
            app.GridLayout3.Layout.Column = 1;

            % Create LastButton
            app.LastButton = uibutton(app.GridLayout3, 'push');
            app.LastButton.ButtonPushedFcn = createCallbackFcn(app, @LastButtonPushed, true);
            app.LastButton.Enable = 'off';
            app.LastButton.Layout.Row = 1;
            app.LastButton.Layout.Column = 1;
            app.LastButton.Text = 'Last';

            % Create NextButton
            app.NextButton = uibutton(app.GridLayout3, 'push');
            app.NextButton.ButtonPushedFcn = createCallbackFcn(app, @NextButtonPushed, true);
            app.NextButton.Enable = 'off';
            app.NextButton.Layout.Row = 1;
            app.NextButton.Layout.Column = 3;
            app.NextButton.Text = 'Next';

            % Create EventsSlider
            app.EventsSlider = uislider(app.GridLayout3);
            app.EventsSlider.Limits = [1 100];
            app.EventsSlider.ValueChangedFcn = createCallbackFcn(app, @EventsSliderValueChanged, true);
            app.EventsSlider.Enable = 'off';
            app.EventsSlider.Layout.Row = 1;
            app.EventsSlider.Layout.Column = 2;
            app.EventsSlider.Value = 1;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout2);
            app.GridLayout4.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout4.RowHeight = {'1x'};
            app.GridLayout4.Layout.Row = 1;
            app.GridLayout4.Layout.Column = 1;

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GridLayout4);
            app.GridLayout5.ColumnWidth = {'1x'};
            app.GridLayout5.Padding = [0 0 0 0];
            app.GridLayout5.Layout.Row = 1;
            app.GridLayout5.Layout.Column = 1;

            % Create TruePositiveRippleLabel
            app.TruePositiveRippleLabel = uilabel(app.GridLayout5);
            app.TruePositiveRippleLabel.HorizontalAlignment = 'center';
            app.TruePositiveRippleLabel.VerticalAlignment = 'top';
            app.TruePositiveRippleLabel.FontWeight = 'bold';
            app.TruePositiveRippleLabel.Layout.Row = 2;
            app.TruePositiveRippleLabel.Layout.Column = 1;
            app.TruePositiveRippleLabel.Text = 'True Positive Ripple';

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout5);
            app.GridLayout7.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout7.RowHeight = {'1x'};
            app.GridLayout7.Padding = [0 0 0 0];
            app.GridLayout7.Layout.Row = 1;
            app.GridLayout7.Layout.Column = 1;

            % Create RippleCheckBox
            app.RippleCheckBox = uicheckbox(app.GridLayout7);
            app.RippleCheckBox.ValueChangedFcn = createCallbackFcn(app, @RippleCheckBoxValueChanged, true);
            app.RippleCheckBox.Enable = 'off';
            app.RippleCheckBox.Text = '';
            app.RippleCheckBox.Layout.Row = 1;
            app.RippleCheckBox.Layout.Column = 4;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout4);
            app.GridLayout6.ColumnWidth = {'1x'};
            app.GridLayout6.Padding = [0 0 0 0];
            app.GridLayout6.Layout.Row = 1;
            app.GridLayout6.Layout.Column = 3;

            % Create TruePositiveSharpWaveLabel
            app.TruePositiveSharpWaveLabel = uilabel(app.GridLayout6);
            app.TruePositiveSharpWaveLabel.HorizontalAlignment = 'center';
            app.TruePositiveSharpWaveLabel.VerticalAlignment = 'top';
            app.TruePositiveSharpWaveLabel.FontWeight = 'bold';
            app.TruePositiveSharpWaveLabel.Layout.Row = 2;
            app.TruePositiveSharpWaveLabel.Layout.Column = 1;
            app.TruePositiveSharpWaveLabel.Text = 'True Positive Sharp Wave';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.GridLayout6);
            app.GridLayout8.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout8.RowHeight = {'1x'};
            app.GridLayout8.Padding = [0 0 0 0];
            app.GridLayout8.Layout.Row = 1;
            app.GridLayout8.Layout.Column = 1;

            % Create SharpwaveCheckBox
            app.SharpwaveCheckBox = uicheckbox(app.GridLayout8);
            app.SharpwaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @SharpwaveCheckBoxValueChanged, true);
            app.SharpwaveCheckBox.Enable = 'off';
            app.SharpwaveCheckBox.Text = '';
            app.SharpwaveCheckBox.Layout.Row = 1;
            app.SharpwaveCheckBox.Layout.Column = 4;

            % Create SaveTruePositiveEventsButton
            app.SaveTruePositiveEventsButton = uibutton(app.GridLayout4, 'push');
            app.SaveTruePositiveEventsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveTruePositiveEventsButtonPushed, true);
            app.SaveTruePositiveEventsButton.WordWrap = 'on';
            app.SaveTruePositiveEventsButton.Enable = 'off';
            app.SaveTruePositiveEventsButton.Layout.Row = 1;
            app.SaveTruePositiveEventsButton.Layout.Column = 2;
            app.SaveTruePositiveEventsButton.Text = 'Save True Positive Events';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ThresholdReviewer_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
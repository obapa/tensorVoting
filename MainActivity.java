package com.example.tensor;

import android.content.Intent;
import android.content.res.Resources;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import static android.provider.AlarmClock.EXTRA_MESSAGE;

public class MainActivity extends AppCompatActivity {
    public static final String id = "com.example.tensor.MESSAGE";
    public static final long tId = 0;
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
    }

    public void startTensor(View view) {

        Intent intent = new Intent(this, calcTensor.class);
        intent.putExtra(id, view.getId());
        startActivity(intent);
/*
        }
        Intent intent = new Intent(this, Tensor.class);
        Button read_dot_edge_file = (Button)findViewById(R.id.TrueLoc1tion);
        //
        Resources res = getResources();
        String[] edge = res.getStringArray(R.array.TrueLocations);
        //
        String message = read_dot_edge_file.getText().toString();
        intent.putExtra(EXTRA_MESSAGE, message);
        startActivity(intent);
*/
    }
}


